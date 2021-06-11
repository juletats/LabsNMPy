#Решение краевых задач для нелинейных дифференциальных уравнений методом конечных разностей.

import numpy as np
from math import pi,cos,sqrt,e
import matplotlib.pyplot as plt

def progon(a,b):#метод прогонки

    # a - исходная трёхдиагоналка
    # b - результат

    n = len(b)-1

    P = [0] * (n+1)
    x = [0] * (n+1)
    Q = [0] * (n+1)

    P[0] = -a[0][1] / a[0][0]
    Q[0] = b[0] / a[0][0]

    for i in range(1,n):
        z = a[i][i] + a[i][i-1] * P[i-1]
        P[i] = -a[i][i+1] / z
        Q[i] = (b[i] - a[i][i-1] * Q[i-1]) / z

    x[n] = (b[n] - a[n][n-1] * Q[n-1]) / (a[n][n] + a[n][n-1] * P[n-1])

    for i in reversed(range(n)):
        x[i] = P[i] * x[i+1] + Q[i]
     
    return x


def f(x):
    return 1-(x*x+1)*(x*x+1)

def KonRaz(X,h,eps = 0.001):
    n = X.shape[0]#количеcтво строк X
    y = [-3 * i for i in X]#пусть начальное приближение равняется 1/3*x
    A = np.zeros((n-2,n-2),dtype = float)#нулевая матрица
    B = np.zeros(n-2,dtype = float)#нулевой вектор
    answer = []#будущие точные значения
    #B[0] = 2*X[0]* (1 + X[0]**3 /3)#то что за знаком равенства
    B[0] = 16*X[0]-4
    #A[0][0] = -2/(h*h) + 2*X[0] * (1 + 0.5*X[0])
    A[0][0] = -2/(h*h)-8
    A[0][1] =  1/(h*h)+y[0]/(2*h*(X[0]*X[0]+1)*(X[0]+2))
    for i in range(1,n-3):
        A[i][i-1] =  1/(h*h)-y[i]/(2*h*(X[i]*X[i]+1)*(X[i]+2))
        A[i][i] = -2/(h*h)-8
        A[i][i+1] =  1/(h*h)+y[i]/(2*h*(X[i]*X[i]+1)*(X[i]+2))
        B[i] = 16*X[i]-4
    A[-1][-2] = 1/(h*h)+y[-2]/(2*h*(X[-2]*X[-2]+1)*(X[-2]+2))
    A[-1][-1] = -2/(h*h)-8
    B[-1] = 16*X[-2]-4+3*1/(h*h)+y[-2]/(2*h*(X[-2]*X[-2]+1)*(X[-2]+2))
    #получили систему с трехдиагональной матрицей, которую решаем методом прогонки
    C=progon(A,B)#решаем методом прогонки
    y_new = [i for i in C]#заполняем новый массив решениями системы
    y_new.insert(0,0)#пусть y0=0
    y_new.append(-3)#пусть yn=1/3
    k = 1

    while True and k < 100:
        y_old = [i for i in y_new]

        A = np.zeros((n-2,n-2),dtype = float)
        B = np.zeros(n-2,dtype = float)
        B[0] =16*X[0]-4
        A[0][0] = -2/(h*h)-8
        A[0][1] = 1/(h*h)+y_old[0]/(2*h*(X[0]*X[0]+1)*(X[0]+2))
        for i in range(1,n-3):
            A[i][i-1] = 1/(h*h)-y_old[i]/(2*h*(X[i]*X[i]+1)*(X[i]+2))
            A[i][i] = -2/(h*h)-8
            A[i][i+1] =  1/(h*h)+y_old[i]/(2*h*(X[i]*X[i]+1)*(X[i]+2))
            B[i] = 16*X[i]-4
        A[-1][-2] =  1/(h*h)+y_old[-2]/(2*h*(X[-2]*X[-2]+1)*(X[-2]+2))
        A[-1][-1] = -2/(h*h)-8
        B[-1] = 16*X[-2]-4+3*1/(h*h)+y_old[-2]/(2*h*(X[-2]*X[-2]+1)*(X[-2]+2))
        C=progon(A,B)#решаем методом прогонки
        
        y_new = [i for i in C]
        y_new.insert(0,0)
        y_new.append(-3)

        e = [(y_old[i] - y_new[i])**2 for i in range(len(y_new))]#проверяем условие окончание итерац процесса
        if sqrt(sum(e)) < eps:
            break
        k += 1
    
    answer = [f(i) for i in X]#точное значение
    eps = [abs(y_new[i] - answer[i]) for i in range(X.shape[0])]#разность между полученным и точным
    #print("Eps = {:6f} \n Количество итераций - {}".format(sqrt(sum(e)),k ))
    print("Количество итераций - {}".format(k))
    draw_k(X,y_new,answer,eps)



def draw_k(X,Y,A,E):
    print('____________________________________________________')
    print('k |   x   |  Точное   |  Полученное | Погрешность   |')
    print('__|_______|___________|_____________|_______________|')
    for i in range(X.shape[0]):
        print("{} | {:.4f}|{:.4f}     | {:.4f}      | {:.4f}        |"\
        .format(i,X[i],Y[i],A[i],E[i]))
    print('___|_______|___________|_____________|_______________|')

    plt.title('Сравнение точного и полученного решений')
    plt.plot(X,Y)#график для точного
    plt.plot(X,A,'mD--')#для полученного
    plt.xlabel('x')
    plt.ylabel('y(x)')
    plt.legend(('Точное','Полученное'))
    plt.grid()
    plt.show()

if __name__ == '__main__':
    X_0 = 0.
    X_k = 1.
    h_1 = .1
    X = np.arange(X_0,X_k+h_1, h_1)#последовательность от x0 до xK с шагом h
    A = KonRaz(X,h_1)