import numpy as np

exp = np.e


def NullMassive(n):
    """Нулевой массив величины n"""
    return [0] * n


def f(x, y, z):
    """Сама задача Коши"""
    return 4 * x * z - (4 * x ** 2 - 3) * y + exp ** (x ** 2)


def exect_f(x):
    """Точное решение"""
    return (exp ** x + exp ** (-x) - 1) * exp ** (x ** 2)


def Euler(a, b, h):
    """Метод Эйлера"""
    x = np.arange(a, b + h, h)
    n = len(x)
    y = np.ones(n)
    z = NullMassive(n)
    for i in range(1, n):
        z[i] = z[i - 1] + h * f(x[i - 1], y[i - 1], z[i - 1])
        y[i] = y[i - 1] + h * z[i - 1]
    return y


def Runge_Kutty(a, b, h):
    """Рунге-Кутт"""
    x = np.arange(a, b + h, h)
    n = len(x)
    K = NullMassive(4)
    L = NullMassive(4)
    y = np.ones(n)
    z = NullMassive(n)
    for i in range(1, n):
        for j in range(1, len(K)):
            K[0] = h * z[i - 1]
            L[0] = h * f(x[i - 1],
                         y[i - 1],
                         z[i - 1])
            K[j] = h * (z[i - 1] + L[j - 1] / 2)
            L[j] = h * f(x[i - 1] + h / 2,
                         y[i - 1] + K[j - 1] / 2,
                         z[i - 1] + L[j - 1] / 2)
        deltay = (K[0] + 2 * K[1] + 2 * K[2] + K[3]) / 6
        deltaz = (L[0] + 2 * L[1] + 2 * L[2] + L[3]) / 6
        y[i] = y[i - 1] + deltay
        z[i] = z[i - 1] + deltaz
    return y, z


def Adams(a, b, h):
    """Метод Адамса"""

    x = np.arange(a, b + h, h)
    n = len(x)
    y = NullMassive(n)
    z = NullMassive(n)
    y_start = Runge_Kutty(a, a + 3 * h, h)[0]
    for i in range(len(y_start)):
        y[i] = y_start[i]

    z_start = Runge_Kutty(a, a + 3 * h, h)[1]
    for i in range(len(z_start)):
        z[i] = z_start[i]

    for i in range(4, n):
        z[i] = (z[i - 1] +
                h * (55 * f(x[i - 1], y[i - 1], z[i - 1]) -
                     59 * f(x[i - 2], y[i - 2], z[i - 2]) +
                     37 * f(x[i - 3], y[i - 3], z[i - 3]) -
                     9 * f(x[i - 4], y[i - 4], z[i - 4]))
                / 24)
        y[i] = y[i - 1] + h * z[i - 1]

    return y


def RunRomRich_method(a, b, h):
    """Рунге-Ромберг-Ричардсон. Возвращает 3 массива для трех методов"""
    x = np.arange(a, b + h, h)
    n = len(x)
    Euler1 = NullMassive(n)
    Runge_Kutty1 = NullMassive(n)
    Adams1 = NullMassive(n)
    Euler_norm = Euler(a, b, h)
    Euler_half = Euler(a, b, h / 2)
    Runge_Kutty_norm = Runge_Kutty(a, b, h)[0]
    Runge_Kutty_half = Runge_Kutty(a, b, h / 2)[0]
    Adams_norm = Adams(a, b, h)
    Adams_half = Adams(a, b, h / 2)
    for i in range(n):
        Euler1[i] = Euler_norm[i] + (Euler_half[i * 2] - Euler_norm[i]) / (1 - 0.5 ** 2)
        Runge_Kutty1[i] = Runge_Kutty_norm[i] + (Runge_Kutty_half[i * 2] - Runge_Kutty_norm[i]) / (1 - 0.5 ** 2)
        Adams1[i] = Adams_norm[i] + (Adams_half[i * 2] - Adams_norm[i]) / (1 - 0.5 ** 2)

    return Euler1, Runge_Kutty1, Adams1


#######################################################################################################################

h = 0.1
a = 0;
b = 1

x = np.arange(a, b + h, h)
y = exect_f(x)  # точное значение функции в точках х

n = len(x)

y_Euler = Euler(a, b, h)
y_Runge_Kutt = Runge_Kutty(a, b, h)[0]
y_Adams = Adams(a, b, h)

y_RunRomRich_Euler = RunRomRich_method(a, b, h)[0]
y_RunRomRich_Runge_Kutt = RunRomRich_method(a, b, h)[1]
y_RunRomRich_Adams = RunRomRich_method(a, b, h)[2]

delta_y_RunRomRich_Euler = [];
delta_y_RunRomRich_Runge_Kutt = [];
delta_y_RunRomRich_Adams = []

for i in range(len(y_RunRomRich_Adams)):
    delta_y_RunRomRich_Euler.append(abs(y_RunRomRich_Euler[i] - y[i]))
    delta_y_RunRomRich_Runge_Kutt.append(abs(y_RunRomRich_Runge_Kutt[i] - y[i]))
    delta_y_RunRomRich_Adams.append(abs(y_RunRomRich_Adams[i] - y[i]))

print('\n------------------------------------------------------------------------------------------')
print('|     |            |       |             |       ||  применяя Рунге-Ромберга-Ричардсона  |')
print('|  x  |  y-точное  | Эйлер | Рунге-Кутта | Адамс || ------------------------------------ |')
print('|     |            |       |             |       ||    Эйлер   | Рунге-Кутта |   Адамс   |')
print('------------------------------------------------------------------------------------------')
for i in range(n):
    print(
        ' {0: < 1.1f} |   {1: < 7.3f}  |{2: < 7.3f}|   {3: < 7.3f}   |{4: < 7.3f}||   {5: < 7.3f}  |   {6: < 7.3f}   |  {7: < 7.3f}  |'.format(
            x[i], y[i], y_Euler[i], y_Runge_Kutt[i], y_Adams[i],
            y_RunRomRich_Euler[i], y_RunRomRich_Runge_Kutt[i], y_RunRomRich_Adams[i]))
print('------------------------------------------------------------------------------------------')
print('                   Разница от точного значения   ||---------------------------------------')
for i in range(n):
    print(
        '                                                 ||   {0: < 7.3f}  |   {1: < 7.3f}   |  {2: < 7.3f}  |'.format(
            delta_y_RunRomRich_Euler[i], delta_y_RunRomRich_Runge_Kutt[i], delta_y_RunRomRich_Adams[i]))

print('------------------------------------------------------------------------------------------')