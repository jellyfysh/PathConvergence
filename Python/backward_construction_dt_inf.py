import math, random

def RhoFree(x1, x2, tau): return math.exp( - (x1 - x2) ** 2 / (2.0 * tau))
random.seed('HoltzmanBeauvalletKrauth')

p = 10
beta = 10; x_beta = 1.0; Deltau = 0.1; sigmaT = 5.0; N = int(beta / Deltau)

gamma = {}
for k in range(1, N + 1): 
    tau = k * Deltau 
    gamma[k] = random.gauss(0.0, math.sqrt(Deltau *  (tau - Deltau) / tau)) 
x = {N: x_beta}
for k in range(N, 0, -1):
    tau = k * Deltau
    Dummy = sigmaT ** 2 * (x[k] * beta - x_beta * tau) / (-sigmaT ** 2 * tau +
    beta * (sigmaT ** 2 + tau))
    x_mean = x[k] * (1.0 - Deltau / tau) + Deltau / tau * Dummy
    x[k - 1] = x_mean + gamma[k]
    print(x[k - 1], 'position at', round(tau - Deltau, 2))
