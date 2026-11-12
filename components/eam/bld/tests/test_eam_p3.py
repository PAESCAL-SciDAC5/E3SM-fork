#!/usr/bin/env python3
"""Check EAM P3 configuration without creating, compiling, or running a case.

Run: python3 components/eam/bld/tests/test_eam_p3.py
Only Perl, Python's standard library, and temporary output directories are needed.
The initial-condition/forcing paths are placeholders; no datasets are read.
"""

from pathlib import Path
import re
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET


BLD = Path(__file__).resolve().parents[1]
SCM_NAMELIST = (
    "&camexp ncdata='initial.nc', start_ymd=20000101, scmlat=36.6, "
    "scmlon=262.5, iopfile='forcing.nc', precip_off=.true. /"
)


def run(command, directory):
    result = subprocess.run(
        [str(arg) for arg in command], cwd=str(directory),
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,
    )
    if result.returncode:
        raise AssertionError(result.stdout)
    return result.stdout


def configure(directory, sgs, scheme, grid='ne4np4', extra=()):
    directory.mkdir()
    return run([
        BLD / 'configure', '-ccsm_seq', '-scam', '-dyn', 'se', '-hgrid', grid,
        '-nlev', '80', '-chem', 'none', '-' + sgs + '_sgs',
        '-microphys', scheme, '-cam_bld', directory,
    ] + list(extra), directory)


def namelist(directory, extra=''):
    run([
        BLD / 'build-namelist', '-config', directory / 'config_cache.xml',
        '-csmdata', '/inputdata', '-use_case', 'scam_generic',
        '-namelist', SCM_NAMELIST + extra, '-dir', directory,
    ], directory)
    return (directory / 'atm_in').read_text()


def micro_nl(text):
    return re.search(r'&micro_nl\b(.*?)\n\s*/', text, re.S | re.I).group(1)


class EamP3Configuration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='test-eam-p3-')
        cls.work = Path(cls.temp.name)
        cls.outputs = {}
        for sgs in ('clubb', 'shoc'):
            for scheme in ('p3', 'eam_p3'):
                directory = cls.work / (sgs + '-' + scheme)
                configure(directory, sgs, scheme)
                cls.outputs[sgs, scheme] = namelist(directory)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_sources_and_sgs_selection(self):
        for sgs in ('clubb', 'shoc'):
            for scheme in ('p3', 'eam_p3'):
                with self.subTest(sgs=sgs, scheme=scheme):
                    directory = self.work / (sgs + '-' + scheme)
                    entries = ET.parse(str(directory / 'config_cache.xml')).getroot()
                    config = {e.attrib['id']: e.attrib.get('value', '') for e in entries.findall('entry')}
                    self.assertEqual(config['microphys'], 'p3')
                    self.assertEqual(config['p3_version'], 'eam' if scheme == 'eam_p3' else 'auto')
                    self.assertEqual(config[sgs + '_sgs'], '1')
                    self.assertEqual(config[('shoc' if sgs == 'clubb' else 'clubb') + '_sgs'], '0')
                    source = 'scream' if (sgs, scheme) == ('shoc', 'p3') else 'eam'
                    paths = (directory / 'Filepath').read_text().splitlines()
                    p3_paths = [p for p in paths if '/physics/p3/' in p]
                    self.assertEqual(len(p3_paths), 1)
                    self.assertTrue(p3_paths[0].endswith('/physics/p3/' + source))
                    self.assertRegex(self.outputs[sgs, scheme], r'precip_off\s*=\s*\.true\.')

    def test_clubb_alias_preserves_build_and_namelist(self):
        for name in ('Filepath', 'CIME_cppdefs', 'atm_in', 'drv_flds_in', 'docn_in'):
            self.assertEqual(
                (self.work / 'clubb-p3' / name).read_text(),
                (self.work / 'clubb-eam_p3' / name).read_text(), name,
            )

    def test_shoc_eam_p3_uses_complete_eam_defaults(self):
        eam = micro_nl(self.outputs['shoc', 'eam_p3'])
        self.assertEqual(eam, micro_nl(self.outputs['clubb', 'eam_p3']))
        self.assertIn('4.1.2', eam)
        self.assertIn('/atm/cam/physprops', eam)
        self.assertRegex(eam, r'p3_nc_autocon_expon\s*=\s*-1.10D0')
        self.assertRegex(eam, r'p3_wbf_coeff\s*=\s*1.0')
        self.assertIn('do_cooper_inp3', eam.lower())
        legacy = micro_nl(self.outputs['shoc', 'p3'])
        self.assertIn('4.1.1', legacy)
        self.assertIn('/atm/scream/tables', legacy)
        self.assertNotIn('p3_autocon_coeff', legacy)
        # Selecting EAM P3 retains the SHOC compile-time physics switches.
        self.assertEqual(
            (self.work / 'shoc-p3/CIME_cppdefs').read_text(),
            (self.work / 'shoc-eam_p3/CIME_cppdefs').read_text(),
        )

    def test_grid_specific_defaults_and_user_override(self):
        for sgs in ('clubb', 'shoc'):
            directory = self.work / (sgs + '-ne120')
            configure(directory, sgs, 'eam_p3', grid='ne120np4')
            output = micro_nl(namelist(directory))
            self.assertRegex(output, r'p3_embryonic_rain_size\s*=\s*0.000020D0')
            overridden = micro_nl(namelist(directory, '&micro_nl p3_wbf_coeff=0.25 /'))
            self.assertRegex(overridden, r'p3_wbf_coeff\s*=\s*0.25')

    def test_old_caches_preserve_namelists(self):
        for sgs in ('clubb', 'shoc'):
            directory = self.work / ('old-' + sgs)
            directory.mkdir()
            cache = ET.parse(str(self.work / (sgs + '-p3') / 'config_cache.xml'))
            for entry in list(cache.getroot()):
                if entry.attrib.get('id') == 'p3_version':
                    cache.getroot().remove(entry)
            cache.write(str(directory / 'config_cache.xml'))
            self.assertEqual(namelist(directory), self.outputs[sgs, 'p3'])

    def test_reject_explicit_eam_p3_with_mmf(self):
        with self.assertRaisesRegex(AssertionError, 'does not support MMF'):
            configure(self.work / 'mmf', 'shoc', 'eam_p3', extra=('-use_MMF',))


if __name__ == '__main__':
    unittest.main(verbosity=2)
