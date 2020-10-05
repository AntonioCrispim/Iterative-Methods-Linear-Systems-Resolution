import numpy as np
import math
def lipkin(n,g,h):
    H=np.zeros((n+1,n+1))
    e=np.zeros((n+1,1))
    for i in range(0,n+1):
        for j in range(0,n+1):
            if i==j:
                H[i,j] = -(1+g)/n*(n*(n+2)/4 - (-n/2+i)**2-n/2)-2*h*(-n/2+i)
            if (i+2)<=n:
                H[i,i+2] = -(1-g)/(2*n)*(math.sqrt(n*(n+2)/4-(-n/2+i)*((-n/2+i)+1))*math.sqrt(n*(n+2)/4-((-n/2+i)+1)*((-n/2+i)+2)))
            if (i-2)>=0:
                H[i,i-2] = -(1-g)/(2*n)*(math.sqrt(n*(n+2)/4-(-n/2+i)*((-n/2+i)-1))*math.sqrt(n*(n+2)/4-((-n/2+i)-1)*((-n/2+i)-2)))
    e=np.linalg.eigvals(H)
    return e
