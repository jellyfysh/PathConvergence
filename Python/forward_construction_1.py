import math, random

def RhoFree(x1, x2, tau): return math.exp( - (x1 - x2) ** 2 / (2.0 * tau))
random.seed('HoltzmanBeauvalletKrauth')

beta = 10; x_zero = 1.0; Deltau = 0.1; sigmaR = 5.0; N = int(beta / Deltau)
gamma = {}
for k in range(N): 
    tau = k * Deltau 
    sigma = math.sqrt(Deltau * (1.0 - Deltau * (beta - sigmaR ** 2) / 
        (beta * (beta - tau) + sigmaR ** 2 * tau )))
    gamma[k] = random.gauss(0.0, sigma)

x = {0: x_zero} 
for k in range(N):
    tau = k * Deltau 
    x_mean = x[k] - Deltau * (x[k] * (beta - sigmaR ** 2) + x[0] * sigmaR ** 2) / \
    (beta * (beta - tau) + sigmaR ** 2 * tau)
    x[k + 1] = x_mean + gamma[k]
    print(x[k + 1], 'position at', round(tau + Deltau, 2))
