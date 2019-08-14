#simplex algorithm

import numpy as np

def sa(B, N, b, cb, cn, xb, xn):
    
    sigma = cn - np.dot(np.dot(cb, np.linalg.inv(B)), N)
    sk = max(sigma)
    
    while sk > 0:
        
        for i in range(0, len(sigma)):
            if sigma[i] == sk:
                k = i
                break
            
        theta = 99999
        for i in range(0, len(b)):
             if N[i][k] > 0 and b[i] / N[i][k] < theta:
                    l = i
                    theta = b[i] / N[i][k]
       
        fa = N[l][k] * 1.0
        b[l] = b[l] / fa
        N[l] = N[l] / fa
        B[l] = B[l] / fa
        
        for i in range(0, len(b)):
            if i != l:
                fa = N[i][k] / N[l][k]
                b[i] = b[i] - fa * b[l]
                B[i] = B[i] - fa * B[l]
                N[i] = N[i] - fa * N[l]
        xb[l], xn[k] = [xn[k], xb[l]]
        for i in range(0, len(b)):
            B[i][l], N[i][k] = [N[i][k], B[i][l]]
        cb[l], cn[k] = [cn[k], cb[l]]
        sigma = cn - np.dot(np.dot(cb, np.linalg.inv(B)), N)
        sk = max(sigma)
       
    return np.dot(cb, b)

M = 99999
N = np.array([[1, -2, 1, 0], [-4, 1, 2, -1], [-2, 0, 1, 0]])
B = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
b = np.array([11, 2, 1])
cb = np.array([0, -M, -M])
cn = np.array([3, -1, -1, 0])
xb = np.array([4, 6, 7])
xn = np.array([1, 2, 3, 5])

print(sa(B, N, b, cb, cn, xb, xn))