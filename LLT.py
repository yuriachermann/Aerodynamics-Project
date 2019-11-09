# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate, integrate
import matplotlib.pyplot as plt


class Wing:
    def __init__(self, c, b, alpha, alphaG, alpha_L_0):
        self.c = c
        self.b = b
        self.S = integrate.quad(self.c, -self.b/2, self.b/2)[0]
        self.AR = self.b**2 / self.S
        self.alpha = alpha
        self.alphaG = alphaG
        self.alpha_L_0 = alpha_L_0


def solve_LLT(wing, N):
    theta = np.linspace(0 + np.pi / (20 * N), np.pi - np.pi / (20 * N), N)
    m = np.zeros((N, N), dtype=float)
    ang = np.zeros(N, dtype=float)

    print(wing.AR)

    for i in range(N):
        for j in range(1, N*2, 2):
            y = -wing.b * np.cos(theta[i]) / 2

            if i == 1:
                m[i, round((j-1)/2)] = j**2
            else:
                m[i, round((j-1)/2)] = np.sin(j * theta[i]) * ((2 * wing.b / (np.pi * wing.c(y))) + (j / np.sin(theta[i])))
            ang[round((j-1)/2)] = (wing.alpha - wing.alpha_L_0 + (2 * y / wing.b) * wing.alphaG) * np.pi / 180

    return np.linalg.solve(m, ang)


def gamma_vec(An, b, Uo, N):
    theta = np.linspace(0 + np.pi / (20 * N), np.pi - np.pi / (20 * N), N)
    aux = np.zeros(N, dtype=float)

    for i in range(N):
        for n in range(N):
            aux[i] = aux[i] + An[n] * np.sin((2 * n + 1) * theta[i])

    gamma = 2 * b * Uo * aux

    return gamma


b = 8
cr = 1.2
# b = 18.6
# cr = 4.6
chord = lambda y: cr * np.sqrt(1 - (2 * y / b) ** 2)
# chord = lambda y: 4.6 - ((4.6-1.9) * y / 9.3)
wing_test = Wing(chord, b, 2, 0, -1.213)
#wing_test = Wing(chord, b, 2, 0, 0)
rho = 1.225
Uo_test = 60
N = 500


An = solve_LLT(wing_test, N)

Cl = An[0] * wing_test.AR * np.pi
sum = 0
for i in range(N):
    sum = (i * 2 + 1) * An[i] ** 2 + sum
CDi1 = wing_test.AR * np.pi * sum
CDi2 = Cl ** 2 / (np.pi * wing_test.AR)

Gamma = gamma_vec(An, wing_test.b, Uo_test, N)

theta = np.linspace(0 + np.pi / (20 * N), np.pi - np.pi / (20 * N), N)
span = -wing_test.b * np.cos(theta) / 2

plt.plot(span, Gamma, color='r', linewidth=0.8)
plt.title("Gamma(y)")
plt.ylabel("Gamma  [m²/s]")
plt.xlabel("span  [m]")
plt.show()
plt.savefig('Gamma(y).png')

A1_analytic = wing_test.alpha * (np.pi / 180) / ((2 * wing_test.b / (np.pi * cr)) + 1)
A1_numeric = An[0]

print("Difference between analytical and numerical A1 is ", abs(A1_analytic - A1_numeric))
