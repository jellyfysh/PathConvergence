import math, random

def RhoFree(x1, x2, tau): return math.exp( - (x1 - x2) ** 2 / (2.0 * tau))
random.seed('HoltzmanBeauvalletKrauth')

p = 10
beta = 10; x_beta = 1.0; Deltau = 0.1; sigmaT = 5.0; N = int(beta / Deltau)

gamma = {}; eta = {}
for k in range(1, N + 1):
    tau = k * Deltau
    gamma[k] = random.gauss(0.0, math.sqrt(Deltau *  (tau - Deltau) / tau))
    eta[k] = random.uniform(0.0, 1.0)
Patterns = []
for mu in range(p):
    Patterns.append(random.gauss(0.0, sigmaT))
x = {N: x_beta}
for k in range(N, 0, -1):
    tau = k * Deltau
    Tower = [RhoFree(Patterns[0], x[k], tau) /
        RhoFree(Patterns[0], x_beta, beta)]
    for mu in range(1, p):
        Tower.append(Tower[mu - 1] + RhoFree(Patterns[mu],x[k], tau) /
               RhoFree(Patterns[mu], x_beta, beta))
    for mu in range(p):
        if eta[k] * Tower[-1] < Tower[mu]: break
    x_mean = (1.0 - Deltau / tau) * x[k] + Deltau / tau * Patterns[mu]
    x[k - 1] = x_mean + gamma[k]
    print(x[k - 1], 'position at', round(tau - Deltau, 2))
