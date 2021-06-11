
def RE(element):
    return round(element, 5)


def y(x):
    return 1 / (x ** 2 + 4)



def Xi(x0, xk, h):
    """Cоздание массива иксов"""
    xi = []
    x_tmp = x0
    i = 0
    while x_tmp != xk:
        xi.append((x0 + h * i))
        i += 1
        x_tmp = xi[-1]
    return xi


def Rectangle(xi, h):
    """Метод прямоугольников"""
    F = [0]
    for i in range(len(xi) - 1):
        F.append(F[i] + (y((xi[i] + xi[i + 1]) / 2)) * h)
    return F


def Trapeze(xi, h):
    """Метод трапеций"""
    F = [0]
    for i in range(len(xi) - 1):
        F.append(F[-1] + ((y(xi[i]) + y(xi[i + 1])) * h) / 2)
    return F


def Sympson(xi, h):
    """Метод Симпсона"""
    F = [0]
    summ = y(xi[0]) + y(xi[-1])
    for i in range(1, len(xi) - 1):
        if i % 2 == 0:
            summ += (2 * y(xi[i]))
        else:
            summ += (4 * y(xi[i]))
            F.append(summ * h / 3)

    return (F, summ * h / 3)


def RungeRomber(Fh, Fkh, p):
    """Рунге-Ромбер"""
    return Fh + ((Fh - Fkh) / (2 ** p - 1))


###########################################################

x0 = -2
xk = 2
h1 = 1
h2 = 0.5
answer = 0.7854


xi1 = Xi(x0, xk, h1)
xi2 = Xi(x0, xk, h2)

F_rec1 = Rectangle(xi1, h1)
F_trap1 = Trapeze(xi1, h1)
(F_sym1, summ) = Sympson(xi1, h1)
F_rec2 = Rectangle(xi2, h2)
F_trap2 = Trapeze(xi2, h2)
(F_sym2, summ) = Sympson(xi2, h2)

# ------------------------------ Шаг h1 ------------------------------

print('\nДля шага h = {0}'.format(h1))
print('\ni |    xi    |    yi    | Rectangle|  Trapeze |  Sympson |')
print('----------------------------------------------------------')
for i in range(len(xi1)):
    if i % 2 == 0:
        print('{0} | {1: < 7.5f} | {2: < 7.5f} | {3: < 7.5f} | {4: < 7.5f} | {5: < 7.5f} |'.format(
            i, xi1[i], y(xi1[i]), F_rec1[i], F_trap1[i], F_sym1[i // 2]))
    else:
        print('{0} | {1: < 7.5f} | {2: < 7.5f} | {3: < 7.5f} | {4: < 7.5f} | -------- |'.format(
            i, xi1[i], y(xi1[i]), F_rec1[i], F_trap1[i]))
print('----------------------------------------------------------\n')

# ------------------------------ Шаг h2 ------------------------------
print('\nДля шага h = {0}'.format(h2))
print('\ni |    xi    |    yi    | Rectangle|  Trapeze |  Sympson |')
print('----------------------------------------------------------')
for i in range(len(xi2)):
    if i % 2 == 0:
        print('{0} | {1: < 7.5f} | {2: < 7.5f} | {3: < 7.5f} | {4: < 7.5f} | {5: < 7.5f} |'.format(
            i, xi2[i], y(xi2[i]), F_rec2[i], F_trap2[i], F_sym2[i // 2]))
    else:
        print('{0} | {1: < 7.5f} | {2: < 7.5f} | {3: < 7.5f} | {4: < 7.5f} | -------- |'.format(
            i, xi2[i], y(xi2[i]), F_rec2[i], F_trap2[i]))
print('----------------------------------------------------------\n')

# ------------------------------ Вывод с Рунге-Робером ------------------------------
print('         Точное решение | {0: < 7.5f}'.format(round(answer, 7)))
print('----------------------------------------------------------')
print('                        | Rectangle|  Trapeze |  Sympson |')
print('----------------------------------------------------------')
print('    Уточненные значения | {0: < 7.5f} | {1: < 7.5f} | {2: < 7.5f} |'.format(
    RungeRomber(F_rec2[-1], F_rec1[-1], 2), RungeRomber(F_trap2[-1], F_trap1[-1], 2),
    RungeRomber(F_sym2[-1], F_sym1[-1], 2)))
print('            Погрешность | {0: < 7.5f} | {1: < 7.5f} | {2: < 7.5f} |'.format(
    abs(answer - RungeRomber(F_rec2[-1], F_rec1[-1], 2)),
    abs(answer - RungeRomber(F_trap2[-1], F_trap1[-1], 2)),
    abs(answer - RungeRomber(F_sym2[-1], F_sym1[-1], 2))))
print('----------------------------------------------------------\n')
print()
