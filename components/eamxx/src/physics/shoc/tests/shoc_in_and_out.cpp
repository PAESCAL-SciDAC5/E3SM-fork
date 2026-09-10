// SHOC "in-and-out" standalone driver.
//
// Reproduces the EAM Fortran turb_standalone ("in-and-out") runs (see
// components/eam/src/physics/cam/shoc_intr.F90, l_turb_standalone path) by
// reading the same three ASCII initial-condition files
//   ShocInOut_IC_surface_vars.txt
//   ShocInOut_IC_zi_grid.txt
//   ShocInOut_IC_zt_grid.txt
// and looping shoc_main over a single column, writing a per-substep dump of
// SHOC's prognostics and diagnostics so the C++ and Fortran engines can be
// compared (BFB) and the dz/dt sensitivity of standalone SHOC studied.
//
// Output: two text tables per run (SHOC uses a staggered vertical grid) —
//   <prefix>_<engine>_<mode>_zt.txt : midpoint (nlev) states/variances vs pres
//   <prefix>_<engine>_<mode>_zi.txt : interface (nlev+1) fluxes/covariances vs presi
// All fluxes are written RAW in native SHOC (kinematic) units; the W/m2
// conversion (x rho x cp / x rho x Lv, matching EAM's shoc_intr history fields)
// is done downstream in the plot script. Values are full double precision.
//
// Two engines are available through the shoc_main_wrap harness:
//   default : the EAMxx C++ shoc_main (SHF::shoc_main via shoc_main_f)
//   -f      : the reference Fortran shoc.F90 compiled into EAMxx
// and two substepping modes matching the EAM namelist l_shoc_outer_loop:
//   exp2 (default, --outer-loop) : N calls to shoc_main with nadv=1,
//                                  per-substep text output after each call
//   exp1 (--single-call)         : one call with nadv=N; only the final
//                                  state can be written (the energy fixer
//                                  runs once per shoc_main call, so exp1
//                                  and exp2 differ physically).
//
// Everything SHOC needs is taken from the IC files or set exactly as
// shoc_intr.F90 sets it (wthv_sec=0, wtracer_sfc=0, tke clipped at mintke,
// phis=0, host_dse = cp*T + g*z). Tunings and physics constants are the
// shared EAM/EAMxx defaults (see shoc_f90.cpp::shoc_init and
// shoc_functions_f90.cpp::shoc_main_f).
//
// High-level flow of main() below:
//   [0] with no arguments, print the usage message and exit
//   [1] READ    : parse the three ShocInOut_IC_*.txt files into an InOutIC
//   [2] PREPARE : convert InOutIC into a FortranData "d" that shoc_main reads,
//                 and initialize the chosen engine (C++ or Fortran)
//   [3] RUN     : loop shoc_main over the substeps, evolving "d" in time
//   [4] WRITE   : after each step, append the column state to the zt & zi tables

#include "shoc_main_wrap.hpp"
#include "shoc_f90.hpp"
#include "shoc_functions_f90.hpp"   // shoc_main_runtime_options()
#include "shoc_constants.hpp"
#include "physics_constants.hpp"

#include "share/scream_types.hpp"
#include "share/scream_session.hpp"

#include "ekat/util/ekat_test_utils.hpp"
#include "ekat/ekat_assert.hpp"

#include <array>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

using scream::Real;
using scream::Int;
// Bring in only the specific SHOC names we use (rather than the whole
// scream::shoc namespace), matching the explicit style of the other tests.
using scream::shoc::FortranData;
using scream::shoc::shoc_init;
using scream::shoc::shoc_main;

// Raw container for the values read out of the three IC text files, kept in
// the files' own "surface-first" order. This is a plain staging struct; the
// flip to SHOC's top-first layout happens later, in make_fortran_data.
struct InOutIC {
  int nz = 0;                                   // number of model levels
  Real dz = 0, zsurf = 0;                       // level thickness, surface height
  Real wqw_sfc = 0, wthl_sfc = 0, uw_sfc = 0, vw_sfc = 0;  // surface fluxes
  // zi-grid profiles, file order (index 0 = surface), size nz+1
  std::vector<Real> zi, presi;
  // zt-grid profiles, file order (index 0 = lowest level), size nz
  std::vector<Real> zt;
  // The 14 zt-grid field blocks, in file order. The enum both names each
  // block and (via NUM_ZT_FIELDS) counts them, so f[THL] is the theta_l block.
  enum ZtField { U=0, V, WM_ZT, TKE, THV, THL, QT, QC, PMID, PDEL,
                 INV_EXNER, TK, TKH, CLOUD_FRAC, NUM_ZT_FIELDS };
  std::array<std::vector<Real>, NUM_ZT_FIELDS> f;
};

// A line counts as a comment (to be skipped) if it is blank or its first
// non-whitespace character is '#'.
bool is_comment (const std::string& line) {
  const auto pos = line.find_first_not_of(" \t\r");
  return pos == std::string::npos || line[pos] == '#';
}

// Read the surface/meta file: 2 comment lines, then one value per line:
// nz, dz, zsurf, wqt_s, wthl_s, t13_s, t23_s (same order shoc_intr reads).
void read_surface_vars (const std::string& fname, InOutIC& ic) {
  std::ifstream in(fname);
  EKAT_REQUIRE_MSG(in, "shoc_in_and_out: cannot open " + fname);
  std::string line;
  std::vector<Real> vals;
  while (std::getline(in, line)) {
    if (is_comment(line)) continue;
    // Parse through a stringstream (with a fail check) so a malformed line
    // gives the same clear message as read_zi_grid/read_zt_grid, instead of
    // throwing a raw std::stod exception.
    std::istringstream ss(line);
    Real val;
    ss >> val;
    EKAT_REQUIRE_MSG(!ss.fail(), "shoc_in_and_out: bad value in " + fname);
    vals.push_back(val);
  }
  EKAT_REQUIRE_MSG(vals.size() == 7,
    "shoc_in_and_out: expected 7 values in " + fname);
  ic.nz       = static_cast<int>(vals[0]);
  ic.dz       = vals[1];
  ic.zsurf    = vals[2];
  ic.wqw_sfc  = vals[3];  // wqt_s  -> wprtp_sfc  [kg/kg m/s]
  ic.wthl_sfc = vals[4];  // wthl_s -> wpthlp_sfc [K m/s]
  ic.uw_sfc   = vals[5];  // t13_s  -> upwp_sfc   [m2/s2]
  ic.vw_sfc   = vals[6];  // t23_s  -> vpwp_sfc   [m2/s2]
}

// Read the zi-grid file: comment lines, then nz+1 rows of "k zi presi",
// surface first (k = 0 at the surface, matching the C++ IC generator).
void read_zi_grid (const std::string& fname, InOutIC& ic) {
  std::ifstream in(fname);
  EKAT_REQUIRE_MSG(in, "shoc_in_and_out: cannot open " + fname);
  std::string line;
  while (std::getline(in, line)) {
    if (is_comment(line)) continue;
    std::istringstream ss(line);
    int k; Real z, p;
    ss >> k >> z >> p;
    EKAT_REQUIRE_MSG(!ss.fail(), "shoc_in_and_out: bad row in " + fname);
    // Rows must arrive in order (k = 0, 1, 2, ...); the running vector size is
    // the index we expect next.
    EKAT_REQUIRE_MSG(k == (int)ic.zi.size(),
      "shoc_in_and_out: non-sequential k index in " + fname);
    ic.zi.push_back(z);
    ic.presi.push_back(p);
  }
  EKAT_REQUIRE_MSG((int)ic.zi.size() == ic.nz+1,
    "shoc_in_and_out: expected nz+1 rows in " + fname);
}

// Read the zt-grid file: a preamble of comment lines, then 14 blocks, each
// introduced by a comment line and holding nz rows of "k zt value", surface
// first. Block order: u, v, wm_zt, tke, thv, thl, qt, qc, pmid, pdel,
// inv_exner, tk, tkh, cloud_frac.
void read_zt_grid (const std::string& fname, InOutIC& ic) {
  std::ifstream in(fname);
  EKAT_REQUIRE_MSG(in, "shoc_in_and_out: cannot open " + fname);
  std::string line;
  // Accumulate rows block by block. Each row is stored as {k, zt, value}.
  std::vector<std::vector<std::array<Real,3>>> blocks;
  while (std::getline(in, line)) {
    if (is_comment(line)) {
      // A comment marks the start of a new block, but only open one if the
      // current block already has data; this collapses runs of consecutive
      // comment lines (the preamble and the per-block headers) into a single
      // boundary.
      if (!blocks.empty() && blocks.back().empty()) continue;
      blocks.emplace_back();
      continue;
    }
    EKAT_REQUIRE_MSG(!blocks.empty(), "shoc_in_and_out: data before header in " + fname);
    std::istringstream ss(line);
    std::array<Real,3> row;
    int k; ss >> k >> row[1] >> row[2];
    row[0] = k;
    EKAT_REQUIRE_MSG(!ss.fail(), "shoc_in_and_out: bad row in " + fname);
    blocks.back().push_back(row);
  }
  // Drop any empty blocks left behind by the preamble comments.
  std::vector<std::vector<std::array<Real,3>>> data;
  for (auto& b : blocks) if (!b.empty()) data.push_back(std::move(b));
  EKAT_REQUIRE_MSG((int)data.size() == InOutIC::NUM_ZT_FIELDS,
    "shoc_in_and_out: expected 14 field blocks in " + fname + ", got " +
    std::to_string(data.size()));
  // Copy each block into ic.f[b], validating that every block shares the same
  // in-order zt column (block 0 defines zt; the rest must match it).
  for (int b = 0; b < InOutIC::NUM_ZT_FIELDS; ++b) {
    EKAT_REQUIRE_MSG((int)data[b].size() == ic.nz,
      "shoc_in_and_out: block " + std::to_string(b) + " in " + fname +
      " has wrong number of rows");
    ic.f[b].resize(ic.nz);
    for (int r = 0; r < ic.nz; ++r) {
      EKAT_REQUIRE_MSG((int)data[b][r][0] == r,
        "shoc_in_and_out: non-sequential k index in " + fname);
      if (b == 0) ic.zt.push_back(data[b][r][1]);
      else EKAT_REQUIRE_MSG(data[b][r][1] == ic.zt[r],
        "shoc_in_and_out: zt mismatch between blocks in " + fname);
      ic.f[b][r] = data[b][r][2];
    }
  }
}

// Populate FortranData from the IC files. The files are surface-first;
// FortranData (like shoc_main in both languages) is top-first, so rows are
// flipped: file row r -> level index nz-1-r (zt) or nz-r (zi).
FortranData::Ptr make_fortran_data (const InOutIC& ic, Real dx, Real dy) {
  using C  = scream::physics::Constants<Real>;
  using SC = scream::shoc::Constants<Real>;

  // Single column, one dummy tracer; nlevi = nlev + 1 interfaces.
  const Int shcol = 1, nlev = ic.nz, nlevi = ic.nz+1, num_qtracers = 1;
  auto dp = std::make_shared<FortranData>(shcol, nlev, nlevi, num_qtracers);
  auto& d = *dp;

  // Column-only (scalar) inputs.
  d.host_dx(0) = dx;
  d.host_dy(0) = dy;
  d.wthl_sfc(0) = ic.wthl_sfc;
  d.wqw_sfc(0)  = ic.wqw_sfc;
  d.uw_sfc(0)   = ic.uw_sfc;
  d.vw_sfc(0)   = ic.vw_sfc;
  d.phis(0)     = 0;                      // flat ocean surface, as in the SCM runs
  d.wtracer_sfc(0,0) = 0;                 // shoc_intr sets tracer fluxes to zero

  // Interface-grid (zi) profiles, flipped to top-first.
  for (Int k = 0; k < nlevi; ++k) {
    const int r = nlevi-1-k;              // flip: k=0 is model top
    d.zi_grid(0,k) = ic.zi[r];
    d.presi(0,k)   = ic.presi[r];
  }

  // Midpoint-grid (zt) profiles, flipped to top-first.
  for (Int k = 0; k < nlev; ++k) {
    const int r = nlev-1-k;               // flip: k=0 is model top
    d.zt_grid(0,k)   = ic.zt[r];
    d.pres(0,k)      = ic.f[InOutIC::PMID][r];
    d.pdel(0,k)      = ic.f[InOutIC::PDEL][r];
    d.inv_exner(0,k) = ic.f[InOutIC::INV_EXNER][r];
    d.thv(0,k)       = ic.f[InOutIC::THV][r];
    d.w_field(0,k)   = ic.f[InOutIC::WM_ZT][r];
    d.u_wind(0,k)    = ic.f[InOutIC::U][r];
    d.v_wind(0,k)    = ic.f[InOutIC::V][r];
    d.thetal(0,k)    = ic.f[InOutIC::THL][r];
    d.qw(0,k)        = ic.f[InOutIC::QT][r];
    d.shoc_ql(0,k)   = ic.f[InOutIC::QC][r];
    d.tk(0,k)        = ic.f[InOutIC::TK][r];
    d.tkh(0,k)       = ic.f[InOutIC::TKH][r];
    d.shoc_cldfrac(0,k) = ic.f[InOutIC::CLOUD_FRAC][r];
    // shoc_intr clips TKE at mintke (tke_tol) before every shoc_main call
    d.tke(0,k)       = std::max<Real>(SC::mintke, ic.f[InOutIC::TKE][r]);
    // Buoyancy flux is reset to zero for a fresh in-and-out run
    // (shoc_intr.F90 turb_standalone path)
    d.wthv_sec(0,k)  = 0;
    d.qtracers(0,k,0) = 0;
    // Dry static energy, as EAM defines state%s with zsurf = phis = 0:
    //   s = cp*T + g*z,  T = thl/inv_exner + (Lv/cp)*ql
    const Real T = d.thetal(0,k)/d.inv_exner(0,k)
                   + (C::LatVap/C::Cpair)*d.shoc_ql(0,k);
    d.host_dse(0,k) = C::Cpair*T + C::gravit*d.zt_grid(0,k);
  }

  return dp;
}

// Per-substep comprehensive SHOC diagnostic dump, split by vertical grid.
// SHOC's fluxes/covariances live on interfaces (nlev+1) while states and
// variances live on midpoints (nlev), so each run writes two tables:
// *_zt.txt (midpoint) and *_zi.txt (interface). All fluxes are written RAW in
// native SHOC (kinematic) units; the W/m2 conversion that matches EAM/Hui is
// applied downstream in the plot script (see README). Values are printed at
// full double precision (%.17e) so the C++ vs Fortran BFB diff is exact and
// Python parses cleanly. Both engines fill every field in `d`
// (shoc_main_wrap.cpp), so these writers are engine-agnostic.

// Midpoint (zt) table: nlev rows, coordinate = pres. Column order matches the
// header written in main().
void write_state_zt (std::FILE* fp, int time_s, const FortranData& d) {
  using C = scream::physics::Constants<Real>;
  for (Int k = 0; k < d.nlev; ++k) {
    const Real qc = d.shoc_ql(0,k);
    const Real qv = d.qw(0,k) - qc;          // vapor = total water - cloud liquid
    // Temperature back out of theta_l; inv_exner is a fixed input, valid each step.
    const Real T  = d.thetal(0,k)/d.inv_exner(0,k) + (C::LatVap/C::Cpair)*qc;
    std::fprintf(fp,
      "%6d %4d %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e"
      " %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
      time_s, (int)k,
      d.zt_grid(0,k), d.pres(0,k), d.u_wind(0,k), d.v_wind(0,k),
      d.thetal(0,k), T, d.qw(0,k), qv, qc,
      d.shoc_cldfrac(0,k), d.tke(0,k), d.shoc_mix(0,k), d.tk(0,k), d.tkh(0,k),
      d.isotropy(0,k), d.brunt(0,k), d.w_sec(0,k), d.shoc_ql2(0,k),
      d.wthv_sec(0,k), d.wqls_sec(0,k));
  }
}

// Interface (zi) table: nlev+1 rows, coordinate = presi. Fluxes/covariances raw.
void write_state_zi (std::FILE* fp, int time_s, const FortranData& d) {
  for (Int k = 0; k < d.nlevi; ++k) {
    std::fprintf(fp,
      "%6d %4d %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
      time_s, (int)k,
      d.zi_grid(0,k), d.presi(0,k),
      d.wthl_sec(0,k), d.wqw_sec(0,k), d.thl_sec(0,k), d.qw_sec(0,k),
      d.qwthl_sec(0,k), d.uw_sec(0,k), d.vw_sec(0,k), d.w3(0,k), d.wtke_sec(0,k));
  }
}

// Guard for options that consume a following value: make sure one exists.
void expect_another_arg (int i, int argc) {
  EKAT_REQUIRE_MSG(i != argc-1, "Expected another cmd-line arg.");
}

} // namespace

// ---------------------------------------------------------------------------
// SHOC tuning parameters.
//  C++ engine: the 12 runtime options SHOC exposes, by their EAMxx names,
//  mapped onto the SHOCRuntime struct shoc_main_f uses. EAM namelist
//  equivalents: shoc_<name> (e.g. shoc_c_diag_3rd_mom, shoc_Ckh).
//  Fortran engine (-f): the same 12 names plus Ckh_s Ckm_s l_inf_const
//  tscale_const Cee_const and the SHOC+MF constants/switches, set by name on
//  the shoc.F90 module variables through the shoc_iso_c bridge (below).
//  EAM namelist equivalents: shoc_<name> / edmf_<name> (do_precip is
//  edmf_do_precipitation).
extern "C" {
  int shoc_set_param_c(const char* name, Real value);   // 0 ok, 1 unknown name
  int shoc_get_param_c(const char* name, Real* value);  // 0 ok, 1 unknown name
}
static const std::vector<std::string>& f90_param_names () {
  static const std::vector<std::string> t = {
    "lambda_low", "lambda_high", "lambda_slope", "lambda_thresh",
    "thl2tune", "qw2tune", "qwthl2tune", "w2tune", "length_fac", "c_diag_3rd_mom",
    "Ckh", "Ckm", "Ckh_s", "Ckm_s", "l_inf_const", "tscale_const", "Cee_const",
    "mf_L0", "mf_ent0", "mf_nup", "mf_a", "mf_b", "mf_c", "mf_a_wcp", "mf_tau_wcp",
    "do_edmf", "do_condensation", "do_precip", "do_mf_diag", "do_wthv_mf",
    "do_dynamic_L", "do_entr_tke", "do_explicit", "do_integral", "do_implicit",
  };
  return t;
}
using RuntimeOpts = scream::shoc::Functions<Real, scream::DefaultDevice>::SHOCRuntime;
static const std::vector<std::pair<std::string, Real RuntimeOpts::*>>& runtime_option_table () {
  static const std::vector<std::pair<std::string, Real RuntimeOpts::*>> t = {
    {"lambda_low",     &RuntimeOpts::lambda_low},
    {"lambda_high",    &RuntimeOpts::lambda_high},
    {"lambda_slope",   &RuntimeOpts::lambda_slope},
    {"lambda_thresh",  &RuntimeOpts::lambda_thresh},
    {"thl2tune",       &RuntimeOpts::thl2tune},
    {"qw2tune",        &RuntimeOpts::qw2tune},
    {"qwthl2tune",     &RuntimeOpts::qwthl2tune},
    {"w2tune",         &RuntimeOpts::w2tune},
    {"length_fac",     &RuntimeOpts::length_fac},
    {"c_diag_3rd_mom", &RuntimeOpts::c_diag_3rd_mom},
    {"Ckh",            &RuntimeOpts::Ckh},
    {"Ckm",            &RuntimeOpts::Ckm},
  };
  return t;
}

// Read "<name> <value>" lines (blank lines and '#' comments ignored; a name may
// appear once) and store them: C++ engine -> the process-wide options that
// shoc_main_f reads; Fortran engine -> the shoc.F90 module variables via
// shoc_set_param_c. Logical MF switches take 0/1 (or .false./.true. spelled as
// 0/1). Unknown names are an error, so a typo cannot silently keep a default.
static void apply_params_file (const std::string& fname, bool use_fortran) {
  std::ifstream in(fname);
  EKAT_REQUIRE_MSG(in, "shoc_in_and_out: cannot open tuning file " + fname);
  auto& opts = scream::shoc::shoc_main_runtime_options();
  std::map<std::string,int> seen;
  std::string line; int lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    const auto hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);
    std::istringstream iss(line);
    std::string name, extra; double value;
    if (!(iss >> name)) continue;                       // blank / comment-only line
    EKAT_REQUIRE_MSG(iss >> value && !(iss >> extra),
                     "shoc_in_and_out: bad line " + std::to_string(lineno) + " in " + fname + ": '" + line + "'");
    bool found = false;
    if (use_fortran) {
      found = (shoc_set_param_c(name.c_str(), static_cast<Real>(value)) == 0);
    } else {
      for (const auto& kv : runtime_option_table()) {
        if (kv.first == name) { opts.*(kv.second) = static_cast<Real>(value); found = true; break; }
      }
    }
    EKAT_REQUIRE_MSG(found, "shoc_in_and_out: unknown SHOC tuning parameter '" + name + "' in " + fname
                            + (use_fortran ? " (Fortran engine)" : " (C++ engine)"));
    EKAT_REQUIRE_MSG(++seen[name] == 1, "shoc_in_and_out: '" + name + "' given twice in " + fname);
  }
}

// Write the effective tuning parameters next to the output tables.
static void write_params_sidecar (const std::string& fname, const std::string& pfile, bool use_fortran) {
  std::FILE* fp = std::fopen(fname.c_str(), "w");
  EKAT_REQUIRE_MSG(fp, "shoc_in_and_out: cannot write " + fname);
  std::fprintf(fp, "# SHOC runtime (tuning) parameters used by this run (%s engine)\n",
               use_fortran ? "Fortran SHOC+MF" : "C++");
  std::fprintf(fp, "# source: %s\n", pfile.empty() ? "EAM defaults (no -p file)" : pfile.c_str());
  if (use_fortran) {
    for (const auto& name : f90_param_names()) {
      Real v = 0;
      EKAT_REQUIRE_MSG(shoc_get_param_c(name.c_str(), &v) == 0, "shoc_in_and_out: shoc_get_param_c failed for " + name);
      std::fprintf(fp, "%-15s %.17g\n", name.c_str(), static_cast<double>(v));
    }
  } else {
    const auto& opts = scream::shoc::shoc_main_runtime_options();
    for (const auto& kv : runtime_option_table())
      std::fprintf(fp, "%-15s %.17g\n", kv.first.c_str(), static_cast<double>(opts.*(kv.second)));
  }
  std::fclose(fp);
}

int main (int argc, char** argv) {
  // [0] No options at all: print the usage message and exit with a non-zero
  //     code (this is a "wrong usage" exit, not a successful run).
  if (argc == 1) {
    std::cout << argv[0] << " [options]\n"
      "SHOC in-and-out driver: single-column SHOC-only run from\n"
      "ShocInOut_IC_*.txt files, mirroring EAM's turb_standalone mode.\n"
      "Options:\n"
      "  -d <dir>      Directory holding ShocInOut_IC_*.txt. Default: '.'\n"
      "  -f            Use the reference Fortran SHOC instead of C++.\n"
      "  -dt <sec>     SHOC timestep. Default: 60.\n"
      "  -s <steps>    Number of substeps. Default: 360 (6 h at dt=60 s).\n"
      "  -x <mode>     exp2 (default): <steps> shoc_main calls with nadv=1,\n"
      "                per-substep output. exp1: one call with nadv=<steps>,\n"
      "                final state only.\n"
      "  -dx <m>       host_dx = host_dy. Default: 100000.\n"
      "  -o <prefix>   Output file prefix. Default: shoc_in_and_out.\n"
      "  -p <file>     SHOC tuning-parameter file: lines of '<name> <value>',\n"
      "                # comments allowed. C++ engine names: lambda_low\n"
      "                lambda_high lambda_slope lambda_thresh thl2tune qw2tune\n"
      "                qwthl2tune w2tune length_fac c_diag_3rd_mom Ckh Ckm.\n"
      "                Fortran engine (-f, SHOC+MF): those plus Ckh_s Ckm_s\n"
      "                l_inf_const tscale_const Cee_const and the MF settings\n"
      "                mf_L0 mf_ent0 mf_nup mf_a mf_b mf_c mf_a_wcp mf_tau_wcp\n"
      "                do_edmf do_condensation do_precip do_mf_diag do_wthv_mf\n"
      "                do_dynamic_L do_entr_tke do_explicit do_integral\n"
      "                do_implicit (logicals as 0/1). Unlisted names keep the\n"
      "                EAM defaults (do_edmf=0: standard SHOC). The effective\n"
      "                values are always written to <prefix>_<engine>_<mode>.params.\n";
    return 1;
  }

  // [1a] Defaults for every option; command-line flags below override them.
  bool use_fortran = false;                 // false => C++ engine, true => Fortran
  int dt = 60, nsteps = 360;                // timestep [s] and number of substeps
  Real dx = 100000.0;                       // horizontal grid size [m]
  std::string icdir = ".", mode = "exp2", prefix = "shoc_in_and_out";
  std::string pfile;                        // -p: tuning-parameter file (optional)

  // [1b] Walk the arguments once. argv_matches does an exact match, so the
  //      order of the branches does not matter. An option that takes a value
  //      consumes the next token via ++i, so it is skipped by the loop and
  //      never re-examined. Anything we do not recognize is a hard error, so
  //      a typo (e.g. "-steps" instead of "-s") fails loudly instead of being
  //      silently dropped and running with the default.
  for (int i = 1; i < argc; ++i) {
    if (ekat::argv_matches(argv[i], "-f", "--fortran")) {
      use_fortran = true;
    } else if (ekat::argv_matches(argv[i], "-d", "--ic-dir")) {
      expect_another_arg(i, argc); icdir = argv[++i];
    } else if (ekat::argv_matches(argv[i], "-dt", "--dt")) {
      expect_another_arg(i, argc); dt = std::atoi(argv[++i]);
    } else if (ekat::argv_matches(argv[i], "-s", "--steps")) {
      expect_another_arg(i, argc); nsteps = std::atoi(argv[++i]);
    } else if (ekat::argv_matches(argv[i], "-x", "--mode")) {
      expect_another_arg(i, argc); mode = argv[++i];
    } else if (ekat::argv_matches(argv[i], "-dx", "--dx")) {
      expect_another_arg(i, argc); dx = std::atof(argv[++i]);
    } else if (ekat::argv_matches(argv[i], "-o", "--out-prefix")) {
      expect_another_arg(i, argc); prefix = argv[++i];
    } else if (ekat::argv_matches(argv[i], "-p", "--params")) {
      expect_another_arg(i, argc); pfile = argv[++i];
    } else {
      EKAT_REQUIRE_MSG(false, std::string("shoc_in_and_out: unknown option '")
                              + argv[i] + "'");
    }
  }
  EKAT_REQUIRE_MSG(mode == "exp1" || mode == "exp2",
                   "shoc_in_and_out: -x must be exp1 or exp2");
  EKAT_REQUIRE_MSG(dt > 0 && nsteps > 0, "shoc_in_and_out: bad dt/steps");

  // Start the EAMxx/Kokkos runtime; the matching finalize is at the end.
  scream::initialize_scream_session(argc, argv); {

    // [1] READ the three IC text files into the staging struct.
    InOutIC ic;
    read_surface_vars(icdir + "/ShocInOut_IC_surface_vars.txt", ic);
    read_zi_grid     (icdir + "/ShocInOut_IC_zi_grid.txt", ic);
    read_zt_grid     (icdir + "/ShocInOut_IC_zt_grid.txt", ic);

    // [2] PREPARE the FortranData shoc_main expects, and set the timestep.
    auto d = make_fortran_data(ic, dx, dx);
    d->dtime = static_cast<Real>(dt);

    // Same shoc_init the BFB tests use: EAM/EAMxx shared constants, npbl=nlev.
    // Passing use_fortran here also tells the harness which engine to call.
    shoc_init(d->nlev, use_fortran);

    // [2b] SHOC tuning parameters: override the EAM defaults from the -p file,
    // if given (C++: shoc_main_f copies the options on every call; Fortran:
    // shoc.F90 module variables). Set once here, before the first shoc_main,
    // they cover the whole run.
    if (!pfile.empty()) apply_params_file(pfile, use_fortran);

    // Open the two output tables (midpoint zt + interface zi), named after the
    // engine and substepping mode, and write their column headers. k = 0 is the
    // model top in both. All fluxes are raw/native units (convert in the plot).
    const std::string engine = use_fortran ? "f90" : "cxx";
    const std::string base = prefix + "_" + engine + "_" + mode;
    const std::string zt_name = base + "_zt.txt";
    const std::string zi_name = base + "_zi.txt";
    std::FILE* fzt = std::fopen(zt_name.c_str(), "w");
    EKAT_REQUIRE_MSG(fzt, "shoc_in_and_out: cannot write " + zt_name);
    std::FILE* fzi = std::fopen(zi_name.c_str(), "w");
    EKAT_REQUIRE_MSG(fzi, "shoc_in_and_out: cannot write " + zi_name);
    std::fprintf(fzt, "time k zt P u v thetal T qw qv qc cldfrac tke shoc_mix"
                      " tk tkh isotropy brunt w_sec ql2 wthv_sec wqls_sec\n");
    std::fprintf(fzi, "time k zi Pi wthl_sec wqw_sec thl_sec qw_sec qwthl_sec"
                      " uw_sec vw_sec w3 wtke_sec\n");

    // Provenance: the tuning parameters this run actually used (defaults or -p).
    write_params_sidecar(base + ".params", pfile, use_fortran);

    std::printf("shoc_in_and_out: engine=%s mode=%s nz=%d dt=%d steps=%d"
                " dx=%g ic=%s\n", engine.c_str(), mode.c_str(), ic.nz, dt,
                nsteps, dx, icdir.c_str());

    // [3]+[4] RUN the engine and WRITE the state.
    if (mode == "exp2") {
      // One shoc_main call per substep (nadv=1). "d" carries the state
      // forward, so each call's output is the next call's input; the column
      // state is written after every step.
      d->nadv = 1;
      for (int istep = 1; istep <= nsteps; ++istep) {
        shoc_main(*d, use_fortran);
        write_state_zt(fzt, istep*dt, *d);
        write_state_zi(fzi, istep*dt, *d);
      }
    } else { // exp1
      // A single shoc_main call that advances nadv=nsteps substeps internally;
      // only the final state is available to write.
      d->nadv = nsteps;
      shoc_main(*d, use_fortran);
      write_state_zt(fzt, nsteps*dt, *d);
      write_state_zi(fzi, nsteps*dt, *d);
    }

    std::fclose(fzt);
    std::fclose(fzi);
    std::printf("shoc_in_and_out: wrote %s and %s (pblh = %.4f m)\n",
                zt_name.c_str(), zi_name.c_str(), d->pblh(0));

  } scream::finalize_scream_session();

  return 0;
}
