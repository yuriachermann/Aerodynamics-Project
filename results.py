import LLT
import numpy as np
import matplotlib.pyplot as plt

# Submarine Spitfire wing values:
b = 11.23
cr = 1.2755
alpha = 2
alpha_g = 0
alpha_0 = -1.2
rho = 1.225
Uo_test = 200
chord = lambda y: cr * np.sqrt(1 - (2 * y / b) ** 2)
wing_test = LLT.Wing(chord, b, alpha, alpha_g, alpha_0)

N = 500  # Number of points analyzed on the wing
n = 15  # Number of points analyzed on the wing

AN = LLT.solve_LLT(wing_test, N)
An = LLT.solve_LLT(wing_test, n)
A1_analytic = (wing_test.alpha - wing_test.alpha_0) * (np.pi / 180) / ((2 * wing_test.b / (np.pi * cr)) + 1)
print("\n\n\nA1 analytically calculated is ", A1_analytic, "\n")
print("A1 numerically calculated with %g points is " % n, An[0], "\n\t( Difference from analytical is ", abs(A1_analytic - An[0]), ")")
print("A1 numerically calculated with %g points is " % N, AN[0], "\n\t( Difference from analytical is ", abs(A1_analytic - AN[0]), ")")


Cl_numeric_N = AN[0] * wing_test.AR * np.pi
Cl_numeric_n = An[0] * wing_test.AR * np.pi
Cl_analytic = A1_analytic * wing_test.AR * np.pi
print("\n\n\nCl analytically calculated is ", Cl_analytic, "\n")
print("Cl numerically calculated with %g points is " % n, Cl_numeric_n, "\n\t( Difference from analytical is ", abs(Cl_analytic - Cl_numeric_n), ")")
print("Cl numerically calculated with %g points is " % N, Cl_numeric_N, "\n\t( Difference from analytical is ", abs(Cl_analytic - Cl_numeric_N), ")")

sum_n = 0
for i in range(n):
    sum_n = (i * 2 + 1) * An[i] ** 2 + sum_n
sum_N = 0
for i in range(N):
    sum_N = (i * 2 + 1) * AN[i] ** 2 + sum_N
CDi_numeric_N = sum_N * wing_test.AR * np.pi
CDi_numeric_n = sum_n * wing_test.AR * np.pi
CDi_analytic = A1_analytic**2 * wing_test.AR * np.pi
print("\n\n\nCDi analytically calculated is ", CDi_analytic, "\n")
print("CDi numerically calculated with %g points is " % n, CDi_numeric_n, "\n\t( Difference from analytical is ", abs(CDi_analytic - CDi_numeric_n), ")")
print("CDi numerically calculated with %g points is " % N, CDi_numeric_N, "\n\t( Difference from analytical is ", abs(CDi_analytic - CDi_numeric_N), ")")

Gamma_n = LLT.gamma(An, wing_test.b, Uo_test, n)
Gamma_N = LLT.gamma(AN, wing_test.b, Uo_test, N)

theta_n = np.linspace(0 + np.pi / (20 * n), np.pi - np.pi / (20 * n), n)
theta_N = np.linspace(0 + np.pi / (20 * N), np.pi - np.pi / (20 * N), N)
span_n = -wing_test.b * np.cos(theta_n) / 2
span_N = -wing_test.b * np.cos(theta_N) / 2

plt.plot(span_n, Gamma_n, color='b', linewidth=0.8)
plt.plot(span_N, Gamma_N, color='r', linewidth=0.8)
plt.title("Gamma(y)")
plt.ylabel("Gamma  [m²/s]")
plt.xlabel("span  [m]")
plt.show()
plt.savefig('Gamma(y).png')
