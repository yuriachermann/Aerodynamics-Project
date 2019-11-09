# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate, integrate
import matplotlib.pyplot as plt


class Wing:
    # This class contains the wing specifications

    def __init__(self, c, b, alpha, alpha_g, alpha_0):
        self.c = c
        self.b = b
        self.S = integrate.quad(self.c, -self.b/2, self.b/2)[0]
        self.AR = self.b**2 / self.S
        self.alpha = alpha
        self.alpha_g = alpha_g
        self.alpha_0 = alpha_0


def solve_LLT(wing, N):
    # solve_LLT computes the lifting line theory parameters
    # The function returns the Fourier series amplitudes
    #
    # ARGUMENTS:
    # wing: contains the wing specs
    # N: Number of points evaluated on the wing
    #
    # EXAMPLE: An = solve_LLT(wing, N)

    theta = np.linspace(0 + np.pi / (2 * N), np.pi - np.pi / (2 * N), N)
    m = np.zeros((N, N), dtype=float)
    ang = np.zeros(N, dtype=float)

    for i in range(N):
        for j in range(1, N*2, 2):
            y = -wing.b * np.cos(theta[i]) / 2

            if i == 1:
                m[i, round((j-1)/2)] = j**2
            else:
                m[i, round((j-1)/2)] = np.sin(j * theta[i]) * ((2 * wing.b/(np.pi * wing.c(y))) + (j/np.sin(theta[i])))
            ang[round((j-1)/2)] = (wing.alpha - wing.alpha_0 + (2 * y / wing.b) * wing.alpha_g) * np.pi / 180

    return np.linalg.solve(m, ang)


def gamma(An, b, Uo, N):
    # gamma_vec computes the wing circulation distribution
    # The function returns the value of gamma in each point of the wing
    #
    # ARGUMENTS:
    # An: Fourier series amplitudes
    # b:  wing span
    # Uo: freestream velocity
    # N:  Number of points evaluated on the wing
    #
    # EXAMPLE: Gamma = gamma(An, b, Uo, N)

    theta = np.linspace(0 + np.pi / (2 * N), np.pi - np.pi / (2 * N), N)
    aux = np.zeros(N, dtype=float)

    for i in range(N):
        for n in range(N):
            aux[i] = aux[i] + An[n] * np.sin((2 * n + 1) * theta[i])

    g = 2 * b * Uo * aux

    return g
