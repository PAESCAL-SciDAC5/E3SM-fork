import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm
import random
import pickle
from joblib import Parallel, delayed
from joblib import parallel_config
import tempfile
import os
from joblib import load, dump
import netCDF4
import sys
import time
import glob

np.set_printoptions(threshold=sys.maxsize)

atol_u = 1.e-10
atol_t = 1.e-10
atol_q = 1.e-10

############################################################################
# This is an offline routine for computing the bulk interface fluxes
# for ocean-atmosphere exchange. This routine is based on the Large
# and Pond formulation and is an offline version of the shr_flux_atmocn
# subroutine in E3SM. We consider the system of equations
#
#		u* = C_DN(u10n) / (1 + C_DN(u10n)/k * (log(z/10) - psim(zeta))) * U
#		u10n = 1 / (1 + C_DN(u10n)/k * (log(z/10) - psim(zeta))) * U
#		theta* = C_HN(zeta) / (1 + C_HN(zeta)/k * (log(z/10) - psih(zeta))) * d(theta)
#		q* = C_EN / (1 + C_EN/k * (log(z/10) - psiq(zeta))) * d(q)
#
# and apply a fixed point iteration to obtain an approximate solution. The 
# iteration takes the form
#
#		x_{n+1} = dt * f(x_{n}) + (1-dt) * x_{n}		(1)
#
# where 0 < dt <= 1 is the relaxation factor. When dt=1, (1) coincides with
# the standard fixed point iteration.
#
# params[in]:
#		dt 	     : relaxation factor (must be between 0 and 1)
#		tFinal   : stopping "time"; time is measured as t_{n+1} = t_{n} + dt
#		myeps    : regularization flag for C_HN:
#					  C_HN = 0.0327, 			if zeta < -myeps
#					   	     0.018 - (0.0327-0.018)/(2*myeps) * (zeta-myeps), 
#											    if -myeps <= zeta <= myeps
#						     0.018,			    if zeta > myeps
#		zetaClip : flag for turning on stability limiter. When True, we take
#					zeta = min{|zeta|, zetaMax} * sgn(zeta)
#		zetaMax  : maximum value for zeta when limiter is on
############################################################################

# parameters for iteration
tFinal = 1000.0
dt = 1.0
myeps = 0.5
zetaClip = True # flag for stability limiter
zetaMax = 10.0  # if zetaClip = True, then this is the max value of zeta
rhnReg = False  # flag for regularization
numOscillating = 0

# physical constants
k = 0.4 	# von karman
g = 9.80616 # grav. acceleration
# zvir = 28.966 / 18.016 - 1.0
zvir = 0.606

# specify save directory and create it if it doesn't exist
saveDir = '/Users/dong9/OneDrive - LLNL/Figures/LargePond_RHSfn'
if not os.path.exists(saveDir):
    os.makedirs(saveDir)
caseDir = '/p/lustre2/dong9/e3sm-scratch/dane/cases/EAMv3_test_baregroundfluxes_sfc_inout_ne30pg2_oECv3'


############################################################################
# manually overwrite first entry of data for testing. the below results are
# direct output from E3SM and can be compared with the results of this
# routine to ensure that the iteration is converging to the correct solution.
#
# Number of iterations performed: 8 
#
# Atmospheric state variables:
# 	Wind speed components:              -5.0096601325718160 and -2.6507466183624988 
# 	Effective wind speed (this is Umag): 5.6677290230426234      
# 	theta_a:    						 298.96449018340951      
# 	theta_s:   						     302.38129799875850      
# 	Density (rbot):   				     1.1805715834673665      
# 	Height:    						     13.297814690034771      
# 	q_a:    							 1.5706012267409919E-002 
# 	q_s:    							 2.4551575304067659E-002 
#
# Output:
# 	ustar:   0.19575574108448721      
#	tstar:  -0.12721966999851686      
#   qstar:  -3.5131779696148842E-004 
#   zeta:   -0.86712068320339730      
#   u10n:    6.0397863094621442, 
#   C_HN:    3.2700000000000000E-002 
#   C_H:     3.7233487182691545E-002
############################################################################
# umag_array[0] = np.sqrt((-5.0096601325718160)**2 + (-2.6507466183624988)**2)
# thetaa_array[0] = 298.96449018340951 
# thetas_array[0] = 302.38129799875850
# rho_array[0] = 1.1805715834673665
# z_array[0] = 13.297814690034771 
# qa_array[0] = 1.5706012267409919e-002

# umag_array[0] = umag_array[25]
# thetaa_array[0] = thetaa_array[25]
# thetas_array[0] = thetas_array[25]
# rho_array[0] = rho_array[25]
# z_array[0] = z_array[25]
# qa_array[0] = qa_array[25]


############################################################################
# function funcZeta: computes the Obukhov stability parameter zeta. Positive
#					 values correspond to stable regime, negative to unstable, and
#					 zero to neutral.
# params[in]:
# 		 ustar 	   : friction velocity
#		 thetastar : temperature scale
#		 qstar     : humidity scale
#		 zetaMax   : maximum value for zeta when limiter is active
#	     z         : atmospheric height
#		 thetaa    : potential temperature of atmosphere
#		 qa        : specific humidity of atmosphere
#		
# params[out]:
#		 zeta   : Obukhov stability parameter
#		 zeta_u : derivative of zeta wrt ustar
#		 zeta_t : derivative of zeta wrt thetastar
#		 zeta_q : derivative of zeta wrt qstar
############################################################################
def funcZeta(ustar, thetastar, qstar, zetaMax, z, thetaa, qa):
    # zeta = k*g*z * (thetastar/thetaa + qstar/(1.0/zvir + qa)) * 1.0/(ustar**2)
    
    # NOTE: this is different than the ocean algorithm
    thva = thetaa * (1.0 + zvir * qa)
    zeta = k*g*z * (thetastar/thva + qstar/(1.0/zvir + qa)) * 1.0/(ustar**2)
    
    zeta_u = -2.0*k*g*z/(ustar**3) * (thetastar/thetaa + qstar/(1.0/zvir + qa))
    zeta_t = k*g*z/(thetaa*ustar**2)
    zeta_q = k*g*z/((1.0/zvir + qa) * ustar**2)

    if zetaClip:
        if (np.abs(zeta) > zetaMax):
            zeta = zetaMax * np.sign(zeta)

            zeta_t = 0.
            zeta_q = 0.
    # zeta = np.minimum(zetaMax, np.abs(zeta)) * np.sign(zeta)
    return zeta, zeta_u, zeta_t, zeta_q


############################################################################
# function funcChi: computes the auxiliary variable chi
#						chi(zeta) = |1 - 16zeta|^(1/4)       (2)
#
# params[in]:
# 		zeta : Obukhov stability parameter
# 
# params[out]:
#		chi 	 : chi from eq (2)
#		chi_zeta : derivative of chi wrt zeta
#
# NOTE: the use of the limiter in the sqrt is done in E3SM; this implementation 
#       is an exact replication of the E3SM subroutine. Its purpose may be 
#       obsolete since zeta is already limited.
############################################################################
def funcChi(zeta):
    # chi = (np.abs(1. - 16.*zeta))**(1/4)
    chi = np.sqrt(max(np.sqrt(np.abs(1.0 - 16.0*zeta)), 1.0))
    chi_zeta = (1.0/4.0) * (np.abs(1. - 16.*zeta))**(-3/4) * np.sign(1. - 16.*zeta) * -16.0

    return chi, chi_zeta


############################################################################
# function funcPsim: computes the stability function psi_m
#						psi_m = log((1 + chi*(2 + chi))(1 + chi^2)/8) - 2arctan(chi) + pi/2,
#														if zeta < 0
#								-5zeta,					if zeta >= 0
#
# params[in]:
# 		zeta : Obukhov stability parameter
#
# params[out]:
#		psim      : value of psi_m at zeta
#		psim_zeta : derivative of psi_m wrt zeta
#
# NOTE: pi/2 in the definition of psi_m is implement as the truncated 1.571
#       in E3SM
############################################################################
def funcPsim(zeta):
    chi, chi_zeta = funcChi(zeta)

   # psimhu(xd) = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) - 2.0_R8*atan(xd) + 1.571_R8
   # psixhu(xd) = 2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)

    if (zeta < 0.):
        # psim = np.log((1.0 + chi*(2.0 + chi))*(1.0 + chi**2)/8) - 2.0*np.arctan(chi) + 1.571
        # psim_zeta = 2.0/(1.0 + chi) * chi_zeta + 2.0*chi/(1.0 + chi**2) * chi_zeta - 2.0/(1.0 + chi**2) * chi_zeta
        
        psim = np.log((1.0 + chi*(2.0 + chi))*(1.0 + chi**2)/8) - 2.0*np.arctan(chi) + 3.14159265358979323846 * 0.5
        psim_zeta = 2.0/(1.0 + chi) * chi_zeta + 2.0*chi/(1.0 + chi**2) * chi_zeta - 2.0/(1.0 + chi**2) * chi_zeta
    else:
        # psim = -5.0*zeta 
        # psim_zeta = -5.0
        
        psim = -(0.7 * zeta + 0.75 * (zeta - 14.3) * np.exp(-0.35 * zeta) + 10.7)
        psim_zeta = -(0.7 + 0.75 * (np.exp(-0.35 * zeta) + (zeta - 14.3) * -0.35 * np.exp(-0.35 * zeta)))

    return psim, psim_zeta

############################################################################
# function funcPsih: computes the stability function psi_h
#						psi_h =  2log((1+chi^2)/2), if zeta < 0
#								 -5zeta, 		    if zeta >= 0
#
# params[in]:
# 		zeta : Obukhov stability parameter
#
# params[out]:
#		psih      : value of psi_h at zeta
#		psih_zeta : derivative of psi_h wrt zeta
#
# NOTE: psi_q = psi_h.
############################################################################
def funcPsih(zeta):
    chi, chi_zeta = funcChi(zeta)

    if (zeta < 0.):
        psih = 2.0*np.log(0.5*(1.0 + chi**2)) 
        psih_zeta = 4.0*chi/(1.0 + chi**2) * chi_zeta
    else:
        # psih = -5.0*zeta
        # psih_zeta = -5.0
        
        psih = -(0.7 * zeta + 0.75 * (zeta - 14.3) * np.exp(-0.35 * zeta) + 10.7)
        psih_zeta = -(0.7 + 0.75 * (np.exp(-0.35 * zeta) + (zeta - 14.3) * -0.35 * np.exp(-0.35 * zeta)))

    return psih, psih_zeta


############################################################################
# function RHS: returns the iteration function f(x):
#
#		f(x) = [RDN(u10n) / (1 + RDN(u10n)/k * (log(z/10) - psim(zeta))) * |U|
#				1 / (1 + RDN(u10n)/k * (log(z/10) - psim(zeta))) * |U|
#				RHN(zeta) / (1 + RHN(zeta)/k * (log(z/10) - psih(zeta))) * d(theta)
#				REN / (1 + REN/k * (log(z/10) - psiq(zeta))) * d(q)]
#
# params[in]:
# 		 ustar 	   : friction velocity
#		 thetastar : temperature scale
#		 qstar     : humidity scale
#		 zetaMax   : maximum value for zeta when limiter is active
#		 Umag      : effective atmospheric wind speed
#	     z         : atmospheric height
#		 thetaa    : potential temperature of atmosphere
#		 thetas    : sea surface temperature
#		 qa        : specific humidity of atmosphere
#		 qs        : specific humidity at sea level
#		 RDN       : exchange coefficient C_HN at previous iteration
#		 rhnReg    : flag for toggling C_HN regularization
#		
# params[out]:
#		 f    : iteration function
#		 RDN  : current value of C_HN
#		 u10n : current 10 m neutral wind speed
#
# NOTE: this function actually just returns the first, third, and fourth 
#		entries of the iteration f. The second entry corresponds to the 
#		neutral 10-meter wind speed u10n and is considered an auxiliary 
#		variable which is not used in the surface flux computations. Thus, 
#		u10n is not an output of this function.
############################################################################
def RHS(ustar, thetastar, qstar, zetaMax, Umag, z, thetaa, thetas, qa, qs, RDN, rhnReg):
    f = np.zeros([3,])

    zeta, _, _, _ = funcZeta(ustar, thetastar, qstar, zetaMax, z, thetaa, qa)
    psim, _ = funcPsim(zeta)
    psih, _ = funcPsih(zeta)

    RHN = RDN 
    REN = RDN
    
    f[0] = RDN / (1.0 + RDN / k * (np.log(z / 10.0) - psim)) * Umag
    f[1] = RHN / (1.0 + RHN / k * (np.log(z / 10.0) - psih)) * (thetaa - thetas)
    f[2] = REN / (1.0 + REN / k * (np.log(z / 10.0) - psih)) * (qa - qs)

    return f


############################################################################
# function eigsRHS: computes eigenvalues of the Jacobian of f(x) at current
#					solution (ustar, thetastar, qstar)
#
# params[in]:
# 		 ustar 	   : friction velocity
#		 thetastar : temperature scale
#		 qstar     : humidity scale
#		 zetaMax   : maximum value for zeta when limiter is active
#		 Umag      : effective atmospheric wind speed
#	     z         : atmospheric height
#		 thetaa    : potential temperature of atmosphere
#		 thetas    : sea surface temperature
#		 qa        : specific humidity of atmosphere
#		 qs        : specific humidity at sea level
#		 RDN       : exchange coefficient C_HN at previous iteration
#		 rhnReg    : flag for toggling C_HN regularization
#		
# params[out]:
#		 eigs  : eigenvalues of Jacobian
#		 normJ : 2-norm of Jacobian
#
# NOTE: this function needs to be modified to account for u10n!
############################################################################
def eigsRHS(ustar, thetastar, qstar, zetaMax, Umag, z, thetaa, thetas, qa, qs, RDN, rhnReg):
    zeta, zeta_u, zeta_t, zeta_q = funcZeta(ustar, thetastar, qstar, zetaMax, z, thetaa, qa)
    psim, psim_zeta = funcPsim(zeta)
    psih, psih_zeta = funcPsih(zeta)

    f = RHS(ustar, thetastar, qstar, zetaMax, Umag, z, thetaa, thetas, qa, qs, RDN, rhnReg)
    J = np.array([[-(f[0]**2)/k * psim_zeta * zeta_u - 0.0, -(f[0]**2)/k * psim_zeta * zeta_t, -(f[0]**2)/k * psim_zeta * zeta_q],
                  [-(f[1]**2)/k * psih_zeta * zeta_u, -(f[1]**2)/k * psih_zeta * zeta_t - 0.0, -(f[1]**2)/k * psih_zeta * zeta_q],
                  [-(f[2]**2)/k * psih_zeta * zeta_u, -(f[2]**2)/k * psih_zeta * zeta_t, -(f[2]**2)/k * psih_zeta * zeta_q - 0.0]])

    if (zeta > myeps):
        RHN = 0.018
        RHN_zeta = 0.0
    elif (zeta > -myeps):
        m = (0.0327 - 0.018)/(-2.0*myeps)
        RHN = 0.018 + m * (zeta - myeps)
        RHN_zeta = m
    else:
        RHN = 0.0327
        RHN_zeta = 0.0


    f2_zeta = ((1.0 + RHN/k * (np.log(z/10.) - psih)) * RHN_zeta - RHN * (RHN_zeta/k * (np.log(z/10.) - psih) * -psih_zeta)) / (1.0 + RHN/k * (np.log(z/10.) - psih))**2

    J[1,0] = f2_zeta * zeta_u
    J[1,1] = f2_zeta * zeta_t - 0.0
    J[1,2] = f2_zeta * zeta_q

    eigs = np.linalg.eigvals(J)
    _, eigvec = np.linalg.eig(J)
    normJ = np.linalg.norm(J, ord=2)

    return eigs, normJ


############################################################################
# function timeStepper: performs the fixed point iteration
#
# params[in]:
# 		 ustar 	   : friction velocity
#		 thetastar : temperature scale
#		 qstar     : humidity scale
#		 zetaMax   : maximum value for zeta when limiter is active
#		 Umag      : effective atmospheric wind speed
#	     z         : atmospheric height
#		 thetaa    : potential temperature of atmosphere
#		 thetas    : sea surface temperature
#		 qa        : specific humidity of atmosphere
#		 qs        : specific humidity at sea level
#		 RDN       : exchange coefficient C_HN at previous iteration
#
# params[out]:
# 		usol     : history for ustar
#		thetasol : " " thetastar
#		qsol     : " " qstar
#	    tsol     : " " "time"
#	    eig1	 : " " eigenvalues
#	    eig2     : " "
#	    eig3     : " "
#	    normJ    : " " norm of Jacobian
#		zeta     : " " stability parameter
#	    res1     : " " residuals
#	    res2     : " "
#		res3     : " "
#	    zeta     : " "
#       i        : number of iterations performed
############################################################################
def timeStepper(ustar, thetastar, qstar, zetaMax, Umag, z, thetaa, thetas, qa, qs, RDN):
    # lists for storing iteration data
    usol = []
    thetasol = []
    qsol = []
    tsol = []
    zeta = []

    eig1 = []
    eig2 = []
    eig3 = []
    normJ = []
    zeta = []

    # append IC to solution array
    usol.append(ustar)
    thetasol.append(thetastar)
    qsol.append(qstar)

    t = 0.0
    tsol.append(t)

    zetai, _, _, _ = funcZeta(ustar, thetastar, qstar, zetaMax, z, thetaa, qa)
    # print("time: ", t, "zeta: ", zetai, "ustar: ", ustar, "thetastar: ", thetastar, "qstar: ,", qstar)

    # eigs, normJi = eigsRHS(ustar, thetastar, qstar, zetaMax, Umag, z, thetaa, thetas, qa, qs, RDN, rhnReg)
    # eig1.append(eigs[0])
    # eig2.append(eigs[1])
    # eig3.append(eigs[2])
    # normJ.append(normJi)
    zeta.append(zetai)

    f = RHS(ustar, thetastar, qstar, zetaMax, Umag, z, thetaa, thetas, qa, qs, RDN, rhnReg)

    # relative residuals
    atol_u = 1.e-10
    atol_t = 1.e-10
    atol_q = 1.e-10

    res1 = np.abs((ustar - f[0]) / (f[0] + atol_u))
    res2 = np.abs((thetastar - f[1]) / (f[1] + atol_t))
    res3 = np.abs((qstar - f[2]) / (f[2] + atol_q))

    i = 0
    while ((res1 > 1.e-13 or res2 > 1.e-13 or res3 > 1.e-13) and 1.e-6 < np.abs(tFinal-t)):
    # while ((res1 > 1.e-13 or res2 > 1.e-13 or res3 > 1.e-13)):    
        zetai, _, _, _ = funcZeta(ustar, thetastar, qstar, zetaMax, z, thetaa, qa)
        zeta.append(zetai)

        f = RHS(ustar, thetastar, qstar, zetaMax, Umag, z, thetaa, thetas, qa, qs, RDN, rhnReg)
        # eigs, normJi = eigsRHS(ustar, thetastar, qstar, zetaMax, Umag, z, thetaa, thetas, qa, qs, deltaq, RDN, False)

        # res1 = np.abs((ustar - f[0]) / f[0])
        # # guard against case when thetastar = 0 (only affects relative residual computation)
        # if (f[1] == 0.0):
        # 	f[1] = 1.e-15
        # res2 = np.abs((thetastar - f[1]) / f[1])
        # res3 = np.abs((qstar - f[2]) / f[2])

        res1 = np.abs((ustar - f[0]) / (f[0] + atol_u))
        res2 = np.abs((thetastar - f[1]) / (f[1] + atol_t))
        res3 = np.abs((qstar - f[2]) / (f[2] + atol_q))

        # perform iteration
        ustar += dt * (f[0] - ustar)
        thetastar += dt * (f[1] - thetastar)
        qstar += dt * (f[2] - qstar)

        # eig1.append(eigs[0])
        # eig2.append(eigs[1])
        # eig3.append(eigs[2])
        # normJ.append(normJi)
        
        t += dt
        tsol.append(t)

        usol.append(ustar)
        thetasol.append(thetastar)
        qsol.append(qstar)

        i += 1
        # if (i % 10 == 0):
        #     print("time: ", t, "zeta: ", zetai, "ustar: ", ustar, "thetastar: ", thetastar, 
        #           "qstar: ,", qstar, ", res1: ", res1, ", res2: ", res2, ", res3: ", res3)

    global numOscillating
    if (res1 > 1.e-3 or res2 > 1.e-3 or res3 > 1.e-3):
        numOscillating += 1
        print("Number of oscillating cells: ", numOscillating)
        # input("Found crashed cell...")

    usol.append(ustar)
    thetasol.append(thetastar)
    qsol.append(qstar)
    zeta.append(zetai)

    return usol, thetasol, qsol, tsol, eig1, eig2, eig3, normJ, zeta, res1, res2, res3, i


############################################################################
# function process: this is a helper function for performing the fixed point 
# 					iteration
#
#
############################################################################
def process(i, n0, forc_hgt_u_array, forc_q_array, qg_array, forc_th_array, t_grnd_array, forc_u_array, 
            forc_v_array, lat_array, lon_array, ugust_array, thm_array, z0mg_array, z0hg_array,
            z0qg_array, forc_rho_array, dataOut):

    forc_hgt_u = forc_hgt_u_array[i]
    forc_q = forc_q_array[i]
    qg = qg_array[i]
    forc_th = forc_th_array[i]
    t_grnd = t_grnd_array[i]
    forc_u = forc_u_array[i]
    forc_v = forc_v_array[i]
    lat = lat_array[i]
    lon = lon_array[i]
    ugust = ugust_array[i]
    thm = thm_array[i]
    z0mg = z0mg_array[i]
    z0hg = z0hg_array[i]
    z0qg = z0qg_array[i]
    forc_rho = forc_rho_array[i]

    # # derived quantities
    # qqq = 11637800.0
    # TTT = 5897.8
    # qsat = qqq * np.exp(-TTT / thetas)
    # qs = qsat / rbot

    # # initial guess for the iteration (NOTE: u10n_0 = Umag)
    # if (CDN > 1e-11):
    #     RDN = np.sqrt(CDN)
    # else:
    #     RDN = 0.4 / np.log(10.0 / 0.0005)
    # RHN = RDN 
    # REN = RDN 
    
    # ustar0 = RDN * Umag
    # thetastar0 = RHN * (thetaa-thetas)
    # qstar0 = REN * (qa - qs)

    # # compute fixed point
    # (usol, thetasol, qsol, tsol, eig1, eig2, eig3, normJ, zeta, res1, res2, res3, 
    # iters) = timeStepper(ustar0, thetastar0, qstar0, zetaMax, Umag, z, thetaa, thetas, 
    #                     qa, qs, RDN)
    # ustar = usol[-1]
    # thetastar = thetasol[-1]
    # qstar = qsol[-1]

    # f = RHS(ustar, thetastar, qstar, zetaMax, Umag, z, thetaa, thetas, qa, qs, RDN, rhnReg)
    # # eigs, _ = eigsRHS(ustar, thetastar, qstar, zetaMax, Umag, z, thetaa, thetas, qa, qs, RDN, rhnReg)

    # # print("Solution ", i, " of ", n0, "- ustar: ", ustar, ", thetastar: ", thetastar, ", qstar: ", qstar, 
    # #     ", res1: ", res1, ", res2: ", res2, ", res3: ", 
    # #     res3, ", Umag: ", Umag_original, ", z: ", z, ", theta_a: ", thetaa, ", theta_s: ", thetas, ", q_s: ", 
    # #     qs, ", q_a: ", qa, ", rbot: ", rbot, ", Num Oscillating: ", numOscillating)


    # store results
    dataOut[i,0] = forc_hgt_u
    dataOut[i,1] = forc_q
    dataOut[i,2] = qg
    dataOut[i,3] = forc_th
    dataOut[i,4] = t_grnd
    dataOut[i,5] = forc_u
    dataOut[i,6] = forc_v 
    dataOut[i,7] = lat
    dataOut[i,8] = lon
    dataOut[i,9] = ugust
    dataOut[i,10] = thm
    dataOut[i,11] = z0mg
    dataOut[i,12] = z0hg
    dataOut[i,13] = z0qg
    dataOut[i,14] = forc_rho
    dataOut[i,15] = 0

    return


def main():

    ############################################################################
    # load data from 1-year E3SM simulation
    # | Variable name          | Description                                                                    |
    # |------------------------|--------------------------------------------------------------------------------|
    # | `cnd01_ICEFRAC_MCTCPL` | Fraction of grid box area covered by sea ice                                   |
    # | `cnd01_LNDFRAC_MCTCPL` | Fraction of grid box area covered by land                                      |
    # | `cnd01_OCNFRAC_MCTCPL` | Fraction of grid box area covered by ocean                                     |
    # | `cnd01_ZBOT_MCTCPL`    | Height above surface at the midpoint of the lowest atmosphere model layer (m)  |
    # | `cnd01_UBOT_MCTCPL`    | Zonal wind at the midpoint of the lowest atmosphere model layer (m/s)          |
    # | `cnd01_VBOT_MCTCPL`    | Meridional wind at the midpoint of the lowest atmosphere model layer (m/s)     |
    # | `cnd01_UGUST_MCTCPL`   | Extra wind speed from gustiness (m/s)                                          |
    # | `cnd01_TBOT_MCTCPL`    | Temperature at the midpoint of the lowest atmosphere model layer (K)           |
    # | `cnd01_THBOT_MCTCPL`   | Potential temperature at the midpoint of the lowest atmosphere model layer (K) |
    # | `cnd01_TS_MCTCPL`      | Surface temperature (K)                                                        |
    # | `cnd01_SST_MCTCPL`     | Sea Surface Temperature (K)                                                    |
    # | `cnd01_QBOT_MCTCPL`    | Specific humidity at the midpoint of the lowest atmosphere model layer (kg/kg) |
    # | `cnd01_SSQ_MCTCPL`     | Specific humidity at the surface (kg/kg)                                       |
    # | `cnd01_SHF_MCTCPL`     | Sensible heat flux  (W/m2)                                                     |
    # | `cnd01_LHF_MCTCPL`     | Latent heat flux   (W/m2)    
    ############################################################################
    months = ['01']
    subsample = 1
    for month in months:
        filenames = glob.glob(caseDir + '/EAMv3_test_baregroundfluxes_sfc_inout_ne30pg2_oECv3.elm.h1.0001-01-01-00000.nc')
        print(filenames)
        # numslices = [30 for i in range(len(filenames))] # number of time slices per file

        # loop over files
        for index, filename in enumerate(sorted(filenames)):
            print("Reading ", filename)
            data = netCDF4.Dataset(filename)

            ncol = data['UBOT'].shape[1]
            nbatch = 40000
            for idx in range(int(ncol / nbatch)):
                left_idx = idx * nbatch 
                right_idx = np.minimum((idx + 1) * nbatch, ncol)

                lat_array = []         # latitudes 
                lon_array = []	       # longitudes
                forc_hgt_u_array = []  # heights of atmospheric layer
                forc_q_array = []	   # specific humidity in atmosphere
                qg_array = []          # specific humidity at sea level (NOT USED)
                forc_th_array = []     # potential temperature of atmosphere
                t_grnd_array = []      # sea surface temperature
                forc_u_array = []      # wind speed, x
                forc_v_array = []      # wind speed, y
                ugust_array = []       # wind gustiness
                thm_array = []   #
                z0mg_array = []
                z0hg_array = []
                z0qg_array = []
                forc_rho_array = []	   # density of air in atmosphere
            
                for timeslice in range(data['time'].shape[0]):
                    # extract relevant data
                    # land_frac = data['cnd01_LNDFRAC_MCTCPL'][timeslice,left_idx:right_idx]
                    # ocn_frac = data['cnd01_OCNFRAC_MCTCPL'][timeslice,left_idx:right_idx]

                    # get indices of ocean cells
                    # idx = np.where(np.abs(ocn_frac - 1.0) < 1.e-15)
                    # idx = idx[0]

                    # print("Number of cells: ", ocn_frac.shape[0], "Number ocean cells: ", len(idx))

                    # filter out ocean data
                    lat = data['latCell'][left_idx:right_idx]
                    lon = data['lonCell'][left_idx:right_idx]
                    forc_hgt_u = data['FORC_HGT_U'][timeslice,left_idx:right_idx]
                    forc_q = data['QBOT'][timeslice,left_idx:right_idx]
                    qg = data['QG'][timeslice,left_idx:right_idx]
                    forc_th = data['THBOT'][timeslice,left_idx:right_idx]
                    t_grnd = data['TG'][timeslice,left_idx:right_idx]
                    forc_u = data['UBOT'][timeslice,left_idx:right_idx]
                    forc_v = data['VBOT'][timeslice,left_idx:right_idx]
                    ugust = data['UGUST'][timeslice,left_idx:right_idx]
                    thm = data['THM'][timeslice,left_idx:right_idx]
                    z0mg = data['Z0MG'][timeslice,left_idx:right_idx]
                    z0hg = data['Z0HG'][timeslice,left_idx:right_idx]
                    z0qg = data['Z0QG'][timeslice,left_idx:right_idx]
                    forc_rho = data['RHOBOT'][timeslice,left_idx:right_idx]
                    
                    print(thm.shape)

                    # location and height
                    lat_array = [*lat_array, *lat[:]]
                    lon_array = [*lon_array, *lon[:]]
                    forc_hgt_u_array = [*forc_hgt_u_array, *forc_hgt_u[:]]
                    forc_q_array = [*forc_q_array, *forc_q[:]]
                    qg_array = [*qg_array, *qg[:]]
                    forc_th_array = [*forc_th_array, *forc_th[:]]
                    t_grnd_array = [*t_grnd_array, *t_grnd[:]]
                    forc_u_array = [*forc_u_array, *forc_u[:]]
                    forc_v_array = [*forc_v_array, *forc_v[:]]
                    ugust_array = [*ugust_array, *ugust[:]]
                    thm_array = [*thm_array, *thm[:]]
                    z0mg_array = [*z0mg_array, *z0mg[:]]
                    z0hg_array = [*z0hg_array, *z0hg[:]]
                    z0qg_array = [*z0qg_array, *z0qg[:]]
                    forc_rho_array = [*forc_rho_array, *forc_rho[:]]

                # global dataOut
                n0 = len(umag_array)
                dataOut = np.zeros([n0, 15])

                print("writing to memmaps")

                temp_folder = tempfile.mkdtemp()
                filename = os.path.join(temp_folder, 'joblib_test.mmap')
                if os.path.exists(filename): os.unlink(filename)
                _ = dump(dataOut, filename)
                large_memmap = load(filename, mmap_mode='r+')

                filename1 = os.path.join(temp_folder, 'joblib_forc_hgt_u.mmap')
                if os.path.exists(filename1): os.unlink(filename1)
                _ = dump(forc_hgt_u_array, filename1)
                forc_hgt_u_memmap = load(filename1, mmap_mode='r+')

                filename2 = os.path.join(temp_folder, 'joblib_forc_q.mmap')
                if os.path.exists(filename2): os.unlink(filename2)
                _ = dump(forc_q_array, filename2)
                forc_q_memmap = load(filename2, mmap_mode='r+')

                filename3 = os.path.join(temp_folder, 'joblib_qg.mmap')
                if os.path.exists(filename3): os.unlink(filename3)
                _ = dump(qg_array, filename3)
                qg_memmap = load(filename3, mmap_mode='r+')

                filename4 = os.path.join(temp_folder, 'joblib_forc_th.mmap')
                if os.path.exists(filename4): os.unlink(filename4)
                _ = dump(forc_th_array, filename4)
                forc_th_memmap = load(filename4, mmap_mode='r+')

                filename5 = os.path.join(temp_folder, 'joblib_t_grnd.mmap')
                if os.path.exists(filename5): os.unlink(filename5)
                _ = dump(t_grnd_array, filename5)
                t_grnd_memmap = load(filename5, mmap_mode='r+')

                filename6 = os.path.join(temp_folder, 'joblib_lat.mmap')
                if os.path.exists(filename6): os.unlink(filename6)
                _ = dump(lat_array, filename6)
                lat_memmap = load(filename6, mmap_mode='r+')

                filename7 = os.path.join(temp_folder, 'joblib_lon.mmap')
                if os.path.exists(filename7): os.unlink(filename7)
                _ = dump(lon_array, filename7)
                lon_memmap = load(filename7, mmap_mode='r+')

                filename8 = os.path.join(temp_folder, 'joblib_forc_u.mmap')
                if os.path.exists(filename8): os.unlink(filename8)
                _ = dump(forc_u_array, filename8)
                forc_u_memmap = load(filename8, mmap_mode='r+')

                filename9 = os.path.join(temp_folder, 'joblib_forc_v.mmap')
                if os.path.exists(filename9): os.unlink(filename9)
                _ = dump(forc_v_array, filename9)
                forc_v_memmap = load(filename9, mmap_mode='r+')
                
                filename10 = os.path.join(temp_folder, 'joblib_ugust.mmap')
                if os.path.exists(filename10): os.unlink(filename10)
                _ = dump(ugust_array, filename10)
                ugust_memmap = load(filename10, mmap_mode='r+')
                
                filename11 = os.path.join(temp_folder, 'joblib_thm.mmap')
                if os.path.exists(filename11): os.unlink(filename11)
                _ = dump(thm_array, filename11)
                thm_memmap = load(filename11, mmap_mode='r+')
                
                filename12 = os.path.join(temp_folder, 'joblib_z0mg.mmap')
                if os.path.exists(filename12): os.unlink(filename12)
                _ = dump(z0mg_array, filename12)
                z0mg_memmap = load(filename12, mmap_mode='r+')
                
                filename13 = os.path.join(temp_folder, 'joblib_z0hg.mmap')
                if os.path.exists(filename13): os.unlink(filename13)
                _ = dump(z0hg_array, filename13)
                z0hg_memmap = load(filename13, mmap_mode='r+')
                
                filename14 = os.path.join(temp_folder, 'joblib_z0qg.mmap')
                if os.path.exists(filename14): os.unlink(filename14)
                _ = dump(z0qg_array, filename14)
                z0qg_memmap = load(filename14, mmap_mode='r+')
                
                filename15 = os.path.join(temp_folder, 'joblib_forc_rho.mmap')
                if os.path.exists(filename15): os.unlink(filename15)
                _ = dump(forc_rho_array, filename15)
                forc_rho_memmap = load(filename15, mmap_mode='r+')

                print("finished writing to memmaps", n0)
                
                # simple parallel implementation; each thread performs the fixed point iteration in
                # the function timeStepper for a single E3SM grid cell in the loaded dataset
                # print(thetaa_memmap)
                with parallel_config(require='sharedmem', n_jobs=4, max_nbytes=None):
                    Parallel()(delayed(process)(i, n0, forc_hgt_u_memmap, forc_q_memmap, qg_memmap, forc_th_memmap, 
                                                t_grnd_memmap, forc_u_memmap, forc_v_memmap, lat_memmap, lon_memmap,
                                                ugust_memmap, thm_memmap, z0mg_memmap, z0hg_memmap, z0qg_memmap, 
                                                forc_rho_memmap, large_memmap) for i in range(n0))
                dataOut = large_memmap

                # save data
                my_dict = {'data': dataOut}
                with open('seaIceOfflineResults_0002-' + month + '_pt' + str(idx) + '.pickle', 'wb') as handle:
                    pickle.dump(my_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    main()

