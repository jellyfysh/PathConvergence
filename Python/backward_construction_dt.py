import math, random

def RhoFree(x1, x2, tau): return math.exp( - (x1 - x2) ** 2 / (2.0 * tau))
random.seed('HoltzmanBeauvalletKrauth')

p = 10
beta = 10; x_beta = 1.0; Deltau = 0.1; sigmaT = 5.0; N = int(beta / Deltau)

gamma = {}
for k in range(1, N + 1): 
    tau = k * Deltau 
    gamma[k] = random.gauss(0.0, math.sqrt(Deltau *  (tau - Deltau) / tau)) 
Patterns = []
for mu in range(p):
    Patterns.append(random.gauss(0.0, sigmaT))
x = {N: x_beta}
for k in range(N, 0, -1):
    tau = k * Deltau
    x_zero_mean = 0.0
    Sum_pi_nu = 0.0
    for nu in range(p):
        pi_nu =  RhoFree(Patterns[nu], x[k], tau) / RhoFree(Patterns[nu], x_beta, beta)
        Sum_pi_nu += pi_nu
        x_zero_mean += pi_nu * Patterns[nu]
    x_mean = x[k] * (1.0 - Deltau / tau) + \
            Deltau / tau * x_zero_mean / Sum_pi_nu
    x[k - 1] = x_mean + gamma[k]
    print(x[k - 1], 'position at', round(tau - Deltau, 2))
