import matplotlib.pyplot as plt
import numpy as np


def NullMatrix(n):
    """Создание нулевой матрицы"""

    nullM = [0] * (n);
    for i in range(n):
        nullM[i] = [0] * (n)
    return nullM


def PrintMatrix(matrix):
    """Печать матрицы"""

    for i in range(len(matrix)):
        print('{0: < 7.3f}  {1: < 7.3f}'.format(round(matrix[i][0], 7), round(matrix[i][1], 7)))


def MultiplyMM(matrix1, matrix2):
    """Умножение матрицы на матрицу"""

    n = len(matrix1);
    result = NullMatrix(n);
    for i in range(n):
        for j in range(n):
            c = 0;
            for k in range(n):
                c += matrix1[i][k] * matrix2[k][j]
            result[i][j] = c
    return result


def MultiplyMV(matrix, vector):
    '''Умножение матрицы на вектор'''

    n = len(vector);
    result = [0] * n;
    for i in range(n):
        for j in range(n):
            result[i] += matrix[i][j] * vector[j];
    return result;


def Swap(matrix, stroke1, stroke2):
    """Смена строк"""

    tmp = 0;
    n = len(matrix)
    for i in range(n):
        tmp = matrix[stroke2][i]
        matrix[stroke2][i] = matrix[stroke1][i]
        matrix[stroke1][i] = tmp;


def SLAE(matrix1, matrix2, vector):
    """Решение СЛАУ"""

    n = len(matrix1)
    vector1 = [0] * n;
    vector2 = [0] * n;

    vector2[0] = vector[0]
    for i in range(1, n):
        summ = 0
        for j in range(0, i):
            summ += matrix1[i][j] * vector2[j]
        vector2[i] = vector[i] - summ
    vector1[n - 1] = vector2[n - 1] / matrix2[n - 1][n - 1]

    for k in range(n - 2, -1, -1):
        summ = 0
        for e in range(k + 1, n):
            summ += matrix2[k][e] * vector1[e];
        vector1[k] = (1 / matrix2[k][k] * (vector2[k] - summ))

    return vector1


def LU(A, b):
    n = len(A)
    M = NullMatrix(n)
    L = NullMatrix(n)
    P = NullMatrix(n)
    P_k = NullMatrix(n)
    for i in range(n):
        for j in range(n):
            if i == j:
                M[i][j] = 1
                L[i][j] = 1
                P[i][j] = 1
                P_k[i][j] = 1

    for i in range(n - 1):
        maximum = -100000


        for l in range(0, n):
            if A[l][i] > maximum:
                maximum = A[l][i]
                index = l;

        Swap(P_k, i, index)

        P = MultiplyMM(P_k, P)
        A = MultiplyMM(P_k, A)

        m = i + 1
        l = i
        for m in range(i + 1, n):
            for l in range(i, n):
                if m == l:
                    M[m][l] = 1
                else:
                    if (m != l) and (l == i):
                        M[m][i] = -A[m][i] / A[i][i]

        A = MultiplyMM(M, A)

        m = i + 1
        for m in range(i + 1, n):
            M[m][i] *= -1

        L = MultiplyMM(P_k, L)
        L = MultiplyMM(L, P_k)

        L = MultiplyMM(L, M)

        M = NullMatrix(n)
        P_k = NullMatrix(n)

        m = 0
        l = 0
        for l in range(0, n):
            for m in range(0, n):
                if l == m:
                    M[l][m] = 1
                    P_k[l][m] = 1

    b = MultiplyMV(P, b)

    x = SLAE(L, A, b)
    return x


x = [-0.7, -0.4, -0.1, 0.2, 0.5, 0.8]
y = [1.3462, 1.9823, 1.671, 1.3694, 1.0472, 0.6435]

SumSquareError1 = 0
SumSquareError2 = 0

sumX = 0;
sumXX = 0;
sumXXX = 0;
sumXXXX = 0
sumY = 0;
sumXY = 0;
sumYXX = 0
for i in range(len(x)):
    sumX += x[i]
    sumXX += x[i] * x[i]
    sumXXX += x[i] * x[i] * x[i]
    sumXXXX += x[i] * x[i] * x[i] * x[i]
    sumY += y[i]
    sumXY += x[i] * y[i]
    sumYXX += x[i] * x[i] * y[i]

N_1 = len(x)  # сразу N+1

matrix = NullMatrix(2)
#  N+1  sumX  = sumY
# sumX  sumXX = sumXY
matrix[0][0] = N_1;
matrix[0][1] = sumX;
matrix[1][0] = sumX;
matrix[1][1] = sumXX;

rightSide = [sumY, sumXY]

A1 = LU(matrix, rightSide)

matrix2 = NullMatrix(3)
matrix2[0][0] = N_1
matrix2[0][1] = sumX
matrix2[0][2] = sumXX
matrix2[1][0] = sumX
matrix2[1][1] = sumXX
matrix2[1][2] = sumXXX
matrix2[2][0] = sumXX
matrix2[2][1] = sumXXX
matrix2[2][2] = sumXXXX

rightSide2 = [sumY, sumXY, sumYXX]

A2 = LU(matrix2, rightSide2)

print('\nПриближающие многочлены')
print('1 степени: F1(x) = {0} + {1}x'.format(round(A1[0], 3), round(A1[1], 3)))
print('2 степени: F2(x) = {0} + {1}x + {2}x^2\n'.format(round(A2[0], 3), round(A2[1], 3), round(A2[2], 3)))

for i in range(len(y)):
    SumSquareError1 += (A1[0] + A1[1] * x[i] - y[i]) ** 2
    SumSquareError2 += (A2[0] + A2[1] * x[i] + A2[2] * x[i] * x[i] - y[i]) ** 2

print('Суммы квадратов ошибок')
print('Ф1 = {0}'.format(round(SumSquareError1, 7)))
print('Ф2 = {0}\n'.format(round(SumSquareError2, 7)))

###########################################################################
X = np.arange(-5, 5, 1)
Y1 = A1[0] + A1[1] * X
Y2 = A2[0] + A2[1] * X + A2[2] * X * X

plt.plot(X, Y1)
plt.plot(X, Y2)
plt.show()