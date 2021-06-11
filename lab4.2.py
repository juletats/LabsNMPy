import numpy as np

def NullMassive(n):
    """Нулевой массив величины n"""
    return [0]*n


def f(x, y, z):
    return -1+2/x

def f_exect(x):
    return np.exp(-x)/x


def p(x):
    return 2/x


def q(x):
    return -1


def Runge_Kutty(a, b, h, y0, z0):
    """Рунге-Кутт"""
    x = np.arange(a, b + h, h)
    K = NullMassive(4)
    L = NullMassive(4)
    y = NullMassive(len(x))
    z = NullMassive(len(x))
    y[0] = y0
    z[0] = z0
    for i in range(1, len(x)):
        for j in range(1, len(K)):
            K[0] = h * z[i - 1]
            L[0] = h * f(x[i - 1], y[i - 1], z[i - 1])
            K[j] = h * (z[i - 1] + L[j - 1] / 2)
            L[j] = h * f(x[i - 1] + h / 2, y[i - 1] + K[j - 1] / 2, z[i - 1] + L[j - 1] / 2)
        deltay = (K[0] + 2 * K[1] + 2 * K[2] + K[3]) / 6
        deltaz = (L[0] + 2 * L[1] + 2 * L[2] + L[3]) / 6
        y[i] = y[i - 1] + deltay
        z[i] = z[i - 1] + deltaz
    return y


def shooting_method(a, b, h, e, y0, y1):
    """Метод стрельбы"""
    nu1 = 1.0
    nu2 = 0.8
    f1 = Runge_Kutty(a, b, h, y0, nu1)[-1] - y1
    f2 = Runge_Kutty(a, b, h, y0, nu2)[-1] - y1

    while(abs(f2) > e):
        nu1, nu2 = nu2, nu2 - f2 * (nu2 - nu1) / (f2 - f1)
        f1, f2 = f2, Runge_Kutty(a, b, h, y0, nu2)[-1] - y1
    return Runge_Kutty(a, b, h, y0, nu2)


def finite_difference_method(a, b, h, alfa, beta, delta, gamma, y0, y1):
    x = np.arange(a, b + h, h)
    N = int((b - a) / h)
    A = []
    B = []
    C = []
    D = []
    A.append(0)
    B.append(-2 + h * h * q(x[1]))
    C.append(1 + p(x[1]) * h / 2)
    D.append(-(1 - (p(x[1]) * h) / 2) * y0)
    for i in range(2, N):
        A.append(1 - p(x[i]) * h / 2)
        B.append(-2 + h * h * q(x[i]))
        C.append(1 + p(x[i]) * h / 2)
        D.append(0)
    A.append(1 - p(x[N - 2]) * h / 2)
    B.append(-2 + h * h * q(x[N - 2]))
    C.append(0)
    D.append(-(1 + (p(x[N - 2]) * h) / 2) * y1)

    P = NullMassive(N)
    Q = NullMassive(N)
    P[0] = (-C[0] / B[0])
    Q[0] = (D[0] / B[0])
    for i in range(1, N):
        P[i] = (-C[i] / (B[i] + A[i] * P[i - 1]))
        Q[i] = ((D[i] - A[i] * Q[i - 1]) / (B[i] + A[i] * P[i - 1]))
    ans = np.zeros(N)
    ans[N - 1] = Q[N - 1]
    for i in range(N - 2, 0, -1):
        ans[i] = P[i] * ans[i + 1] + Q[i]
    ans[0] = y0
    ans = np.append(ans, y1)
    return ans

def RunRomRich_method(a, b, h, e, y0, y1, alfa, beta, gamma, delta):
    x = np.arange(a, b + h, h)
    shooting_method1 = NullMassive(len(x))
    finite_difference_method1 = NullMassive(len(x))
    shooting_method_norm = shooting_method(a, b, h, e, y0, y1)
    shooting_method_half = shooting_method(a, b, h / 2, e, y0, y1)
    finite_difference_method_norm = finite_difference_method(a, b, h, alfa, beta, delta, gamma, y0, y1)
    finite_difference_method_half = finite_difference_method(a, b, h / 2, alfa, beta, delta, gamma, y0, y1)
    for i in range(len(x)):
        shooting_method1[i] = shooting_method_norm[i] + (shooting_method_half[i * 2] - shooting_method_norm[i]) / (1 - 0.5**2)
        finite_difference_method1[i] = finite_difference_method_norm[i] + (finite_difference_method_half[i * 2] - finite_difference_method_norm[i]) / (1 - 0.5**2)
    return shooting_method1, finite_difference_method1
# родные
a = 1; b = 2
h = 0.1
e = 0.001
alfa = 0;  beta = 1;  delta = 0;  gamma = 1
y0 = np.exp(-1)
y1 = np.exp(-2)*0.5

x = np.arange(a, b + h, h)
y = f_exect(x)

n = len(x)

Shooting = shooting_method(a, b, h, e, y0, y1)
Finite = finite_difference_method(a, b, h, alfa, beta, delta, gamma, y0, y1)

Shooting_RunRomRich = RunRomRich_method(a, b, h, e, y0, y1, alfa, beta, delta, gamma)[0]
Finite_RunRomRich = RunRomRich_method(a, b, h, e, y0, y1, alfa, beta, delta, gamma)[1]

delta_shooting = abs(RunRomRich_method(a, b, h, e, y0, y1, alfa, beta, delta, gamma)[0] - y)
delta_Finite = RunRomRich_method(a, b, h, e, y0, y1, alfa, beta, delta, gamma)[1] - y
print('\n-------------------------------------------------------------------------------------')
print('|     |            |          |  Конечно-  ||  применяя Рунге-Ромберга-Ричардсона  |')
print('|  x  |  y-точное  | Стрельба | разностный || ------------------------------------ |')
print('|     |            |          |            ||    Стрельба     | Конечно-разностный |')
print('------------------------------------------------------------------------------------')
for i in range(n):
    print(' {0: < 1.1f} |   {1: < 7.3f}  |  {2: < 7.3f} |   {3: < 7.3f}  ||     {4: < 7.3f}     |       {5: < 7.3f}      |'.format(
        x[i], y[i], Shooting[i], Finite[i], Shooting_RunRomRich[i], Finite_RunRomRich[i]))
print('------------------------------------------------------------------------------------')
print('                   Разница от точного значения   ||---------------------------------')
for i in range(n):
    print('                                                 ||   {0: < 7.3f} |       {1: < 7.3f}      |'.format(
        delta_shooting[i], delta_Finite[i]))

print('------------------------------------------------------------------------------------')