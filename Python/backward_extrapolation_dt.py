import math, random

def RhoFree(x1, x2, tau): return math.exp( -(x1 - x2) ** 2 / (2.0 * tau))
random.seed('HoltzmanBeauvalletKrauth')


ppq = 100 # p plus q
tau_ex = 3
p = ppq // 2
beta = 10; x_beta = 1.0; Deltau = 0.1; sigmaT = 5.0; N = int(beta / Deltau)

gamma = {}
for iter in range(100000):
    Patternsp = []
    Patternsq = []
    for iter in range(p):
        Patternsp.append(random.gauss(0.0, sigmaT))
        Patternsq.append(random.gauss(0.0, sigmaT))
    Patternspq = Patternsp + Patternsq

    pathp = {N: x_beta}; pathq = {N: x_beta}; pathpq = {N: x_beta}; path_inf = {N: x_beta}

    for k in range(N, 0, -1):
        tau = k * Deltau
        gamma[k] = random.gauss(0.0, math.sqrt(Deltau * (tau - Deltau) / tau))
        x_zero_meanp = 0.0; x_zero_meanq = 0.0; x_zero_meanpq = 0.0
        Sum_pi_nup = 0.0; Sum_pi_nuq = 0.0; Sum_pi_nupq = 0.0
        for nu in range(ppq):
            pi_nupq =  RhoFree(Patternspq[nu], pathpq[k], tau) / RhoFree(Patternspq[nu], x_beta, beta)
            Sum_pi_nupq += pi_nupq
            x_zero_meanpq += pi_nupq * Patternspq[nu]
        for nu in range(p):
            pi_nup =  RhoFree(Patternsp[nu], pathp[k], tau) / RhoFree(Patternsp[nu], x_beta, beta)
            pi_nuq =  RhoFree(Patternsq[nu], pathq[k], tau) / RhoFree(Patternsq[nu], x_beta, beta)
            Sum_pi_nup += pi_nup; Sum_pi_nuq += pi_nuq
            x_zero_meanp += pi_nup * Patternsp[nu];  x_zero_meanq += pi_nuq * Patternsq[nu]

        x_tau_minus_deltau_meanp = pathp[k] * (1.0 - Deltau / tau) + Deltau / tau * x_zero_meanp / Sum_pi_nup
        x_tau_minus_deltau_meanq = pathq[k] * (1.0 - Deltau / tau) + Deltau / tau * x_zero_meanq / Sum_pi_nuq
        x_tau_minus_deltau_meanpq = pathpq[k] * (1.0 - Deltau / tau) + Deltau / tau * x_zero_meanpq / Sum_pi_nupq

        Dummy = sigmaT ** 2 * (path_inf[k] * beta - x_beta * tau) /  \
              (- sigmaT ** 2 * tau + beta * (sigmaT ** 2 + tau))
        x_tau_minus_deltau_mean_inf = path_inf[k] * (1 - Deltau / tau) + Deltau / tau * Dummy

        path_inf[k - 1] = x_tau_minus_deltau_mean_inf + gamma[k]
        pathp[k - 1] = x_tau_minus_deltau_meanp + gamma[k]
        pathq[k - 1] = x_tau_minus_deltau_meanq + gamma[k]
        pathpq[k - 1] = x_tau_minus_deltau_meanpq + gamma[k]
    if pathpq[tau_ex] < path_inf[tau_ex]:
        delxp = min(pathp[tau_ex], pathq[tau_ex]) - pathpq[tau_ex]
        delxq = max(pathp[tau_ex], pathq[tau_ex]) - pathpq[tau_ex]
        if -delxp < delxq: print(delxq, -delxp, "o")
        else: print(delxq, -delxp, "x")
