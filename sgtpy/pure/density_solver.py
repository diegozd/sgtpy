import numpy as np
from ..constants import Na
from scipy.optimize import minimize_scalar, brentq


def dPsaft_fun(rho, temp_aux, saft):
    rhomolecular = Na * rho
    global Xass
    da, Xass = saft.d2afcn_aux(rhomolecular, temp_aux, Xass)
    afcn, dafcn, d2afcn = da
    dPsaft = 2 * rhomolecular * dafcn + rhomolecular**2 * d2afcn
    return dPsaft


def Psaft_obj(rho, temp_aux, saft, Pspec):
    rhomolecular = Na * rho
    global Xass
    da, Xass = saft.dafcn_aux(rhomolecular, temp_aux, Xass)
    afcn, dafcn = da
    Psaft = rhomolecular**2 * dafcn / Na
    return Psaft - Pspec


def density_topliss(state, temp_aux, P, Xass0, saft):

    beta = temp_aux[0]
    # lower boundary a zero density
    rho_lb = 1e-5
    dP_lb = Na / beta

    # upper boundary limit at infinity pressure
    etamax = 0.7405
    rho_lim = (6 * etamax) / (saft.ms * np.pi * saft.sigma**3) / Na
    ub_sucess = False
    rho_ub = 0.4 * rho_lim
    it = 0
    P_ub, dP_ub, Xass_ub = saft.dP_drho_aux(rho_ub, temp_aux, Xass0)
    while not ub_sucess and it < 5:
        it += 1
        P_ub, dP_ub, Xass_ub = saft.dP_drho_aux(rho_ub, temp_aux, Xass_ub)
        rho_ub += 0.1 * rho_lim
        ub_sucess = P_ub > P and dP_ub > 0

    # derivative computation  at zero density
    rho_lb1 = 1e-4 * rho_lim
    P_lb1, dP_lb1, Xass_lb = saft.dP_drho_aux(rho_lb1, temp_aux, Xass0)
    d2P_lb1 = (dP_lb1 - dP_lb) / rho_lb1
    if d2P_lb1 > 0:
        flag = 3
    else:
        flag = 1

    global Xass
    Xass = Xass0

    # Stage 1
    bracket = [rho_lb, rho_ub]
    if flag == 1:
        # Found inflexion point
        sol_inf = minimize_scalar(dPsaft_fun, args=(temp_aux, saft),
                                  bounds=bracket, method='Bounded')
        rho_inf = sol_inf.x
        dP_inf = sol_inf.fun
        if dP_inf > 0:
            flag = 3
        else:
            flag = 2

    # Stage 2
    if flag == 2:
        if state == 'L':
            bracket[0] = rho_inf
        elif state == 'V':
            bracket[1] = rho_inf
        rho_ext = brentq(dPsaft_fun, bracket[0], bracket[1],
                         args=(temp_aux, saft))
        P_ext, dP_ext, Xass = saft.dP_drho_aux(rho_ext, temp_aux, Xass)
        if P_ext > P and state == 'V':
            bracket[1] = rho_ext
        elif P_ext < P and state == 'L':
            bracket[0] = rho_ext
        else:
            flag = -1

    if flag == -1:
        rho = np.nan
    else:
        rho = brentq(Psaft_obj, bracket[0], bracket[1],
                     args=(temp_aux, saft, P))

    return rho, Xass


def density_newton(rho0, temp_aux, P, Xass0, saft):

    rho = 1.*rho0
    Psaft, dPsaft, Xass = saft.dP_drho_aux(rho, temp_aux, Xass0)
    FO = Psaft - P
    dFO = dPsaft
    for i in range(15):
        rho -= FO/dFO
        Psaft, dPsaft, Xass = saft.dP_drho_aux(rho, temp_aux, Xass)
        FO = Psaft - P
        dFO = dPsaft
        if np.abs(FO) < 1e-6:
            break
    return rho, Xass
