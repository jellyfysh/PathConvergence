import math, random

def RhoFree(x1, x2, tau): return math.exp( - (x1 - x2) ** 2 / (2.0 * tau))
random.seed('HoltzmanBeauvalletKrauth')

beta = 10; x_zero = 1.0; Deltau = 0.1; sigmaR = 5.0; N = int(beta / Deltau)
gamma = {}
eta = {}
for k in range(N): 
    tau = k * Deltau 
    sigma = 1.0 / math.sqrt(1.0 /(sigmaR ** 2) + 1.0 / (beta - tau) - 1.0 / beta)
    eta[k] = random.gauss(0.0, sigma)
    print(k, (Deltau * (beta - tau - Deltau) / (beta - tau)))
    sigma = math.sqrt(Deltau * (beta - tau - Deltau) / (beta - tau))
    gamma[k] = random.gauss(0.0, sigma)
x = {0: x_zero} 
for k in range(N):
    tau = k * Deltau 
    x_mean = sigmaR ** 2 * ( beta * (x[k] - x[0]) + tau * x[0]) / (beta * (beta
       - tau)  + sigmaR ** 2 * tau)
    xbeta = x_mean + eta[k]
    x_mean = (x[k] * (beta - tau - Deltau) + xbeta * Deltau) / (beta - tau)
    x[k + 1] = x_mean + gamma[k]
    print(x[k + 1], 'position at', round(tau + Deltau, 2))
