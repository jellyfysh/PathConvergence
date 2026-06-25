import math, random

def RhoFree(x1, x2, tau): return math.exp( - (x1 - x2) ** 2 / (2.0 * tau))
random.seed('HoltzmanBeauvalletKrauth')

beta = 10; x_0 = 0.0; x_beta = 1.0; sigmaT = 5.0
x = {0: x_0, beta: x_beta}

for iter in range(16):
    T = sorted(x.keys())
    for tau_minus  in T[:-1]:
        tau_plus = T[T.index(tau_minus) + 1]
        tau = (tau_minus + tau_plus) / 2.0
        x_mean = (x[tau_minus] + x[tau_plus] ) / 2.0
        sigma = math.sqrt(tau_plus - tau_minus) / 2.0
        x[tau] = x_mean + random.gauss(0.0, sigma)
T = sorted(x.keys())
for tau in T:
    print(tau, x[tau])
print(len(T), 'number of points')
