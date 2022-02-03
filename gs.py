import numpy as np

def gaus_seidel(A,b,x0,itMax=100,errMax=0.001):
    n, _  =A.shape
    x=np.zeros(n)
    
    for i in range(itMax):
        for j in range(n):
            x[j]=1/A[j,j] * (b[j]-np.dot(A[j,:j],x[:j])-np.dot(A[j,j+1:],x0[j+1:]))
    
        if np.all(np.abs(x0-x)<errMax):
            print(i)
            break
        x0=x.copy()
    
    return x