import numpy as np
def QR(A):
    a=A.copy()
    Q=np.eye(n)
    for i in range(n-1):
        v=np.zeros(n)
        t=0
        for j in range(i,n):
            t+=(a[j,i])**2
        sign=-1 if a[i,i]<0 else 1
        v[i]=a[i,i]+sign*np.sqrt(t)
        for j in range(i+1,n):
            v[j]=a[j,i]
        p=v.reshape(n,1)*v
        k=v.dot(v.reshape(n,1))
        h=np.eye(n)-2/k[0]*p
        a=h.dot(a)
        Q=Q.dot(h)
    return Q,a
def sum_sqr(a):
    n=a.shape[0]
    s=0
    for i in range(1,n):
        for j in range(0,i):
            s+=a[i,j]**2
    return np.sqrt(s)

array=[[-5, -8, 4],
       [4, 2, 6],
       [-2, 5, 6]]
n=3
A=np.array(array)
print(A)
eps=0.0001
Q, R=QR(A)
print("Проверка QR-алгоритма")
print(A-Q.dot(R))
a=A.copy()
max_it=100
i=0
for i in range(max_it):
    Q,R=QR(a)
    a=R.dot(Q)
    if sum_sqr(a)<eps:
        print(sum_sqr(a))
        break
vec=np.zeros(n, dtype=complex)
i=0
while i<n:
    s=0
    for j in range(i+1,n):
        s+=a[j,i]**2
    if np.sqrt(s)<eps:
        vec[i]=a[i,i]
        i+=1
    else:
        a1=1
        b1=-a[i+1,i+1]-a[i,i]
        c1=a[i,i]*a[i+1,i+1]-a[i,i+1]*a[i+1,i]
        vec[i]=complex(-b1/(2*a1), np.sqrt(abs(b1**2-4*a1*c1))/(2*a1))
        vec[i+1]=complex(-b1/(2*a1), -np.sqrt(abs(b1**2-4*a1*c1))/(2*a1))
        i+=2
for i in range(n):
    print(vec[i])
