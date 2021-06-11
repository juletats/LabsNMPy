def NullMatrix(n):
    """Создание нулевой матрицы"""

    nullM = [0] * (n);
    for i in range(n):
        nullM[i] = [0] * (n)
    return nullM


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
    """Умножение матрицы на вектор"""

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


####################################################################

xi = [-0.4, -0.1, 0.2, 0.5, 0.8]
fi = [1.9823, 1.6710, 1.3694, 1.0472, 0.64350]
X = 0.1

hi = []

# xi = [0,      1,      2,      3,      4]
# fi = [0, 1.8415, 2.9093, 3.1411, 3.2432]
# X = 1.5

n = len(xi)
for i in range(n - 1):
    hi.append(xi[i + 1] - xi[i])

C = NullMatrix(n - 2)
C[0][0] = 2 * (hi[0] + hi[1])
C[0][1] = hi[1]  # i=2
C[1][0] = hi[1]
C[1][1] = 2 * (hi[1] + hi[2])
C[1][2] = hi[2]
C[2][1] = hi[2]
C[2][2] = 2 * (hi[2] + hi[3])

rightSide = []
rightSide.append(3 * (((fi[2] - fi[1]) / hi[1]) - ((fi[1] - fi[0]) / hi[0])))
rightSide.append(3 * (((fi[3] - fi[2]) / hi[2]) - ((fi[2] - fi[1]) / hi[1])))
rightSide.append(3 * (((fi[-1] - fi[-2]) / hi[3]) - ((fi[-2] - fi[-3]) / hi[2])))

ai = [];
bi = [];
di = []
ci = LU(C, rightSide);
ci.insert(0, 0)

for i in range(len(ci) - 1):
    ai.append(fi[i])
    bi.append(((fi[i + 1] - fi[i - 1 + 1]) / hi[i]) - (1 / 3) * hi[i] * (ci[i + 1] + 2 * ci[i]))
    di.append((ci[i + 1] - ci[i]) / (3 * hi[i]))
ai.append(fi[-2])
bi.append((fi[-1] - fi[-2]) / (hi[-1]) - (2 / 3) * hi[-1] * ci[-1])
di.append(-ci[-1] / (3 * hi[-1]))

print('\nКубический сплайн\n')

print('| --- a --- | --- b --- | --- c --- | --- d --- |')
print('| --------------------------------------------- |')
for i in range(len(ai)):
    print('| {0: 9f} | {1: 9f} | {2: 9f} | {3: 9f} |'.format(round(ai[i], 5), round(bi[i], 5), round(ci[i], 5),
                                                             round(di[i], 5)))
print('| --------------------------------------------- |\n')

index = 0
while X > xi[index]:
    index += 1
index -= 1

print('f = {0} + {1}*(x-({2})) + {3}*(x-({4}))^2 + {5}*(x-({6}))^3'.format(round(ai[index], 5), round(bi[index], 5),
                                                                           round(xi[index], 5), round(ci[index], 5),
                                                                           round(xi[index], 5), round(di[index], 5),
                                                                           round(xi[index], 5)))
f = ai[index] + bi[index] * (X - xi[index]) + ci[index] * (X - xi[index]) ** 2 + di[index - 1] * (X - xi[index]) ** 3
print('\nf({0}) = {1}\n'.format(X, round(f, 4)))
