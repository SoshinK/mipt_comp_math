import numpy as np
from scipy.fftpack import dst, idst


n = 100
N = (n - 1) ** 2

def u(x, y):
    return (x-1)*(y-1)*x*y*np.sin(x)

def f(x, y):
    return (x**2 * (y**2-y-2)+x*(-y**2+y+2)-2*(y-1)*y)*np.sin(x) \
           - 2*(2*x-1)*(y-1)*y*np.cos(x)

h = 1. / n

f_h = np.zeros(N)

D = np.zeros((N, N))

for k in range(n - 1):
    for m in range(n - 1):
        f_h[k * (n - 1) + m] = f((k + 1) * h, (m + 1) * h)
        D[k * (n - 1) + m, k * (n - 1) + m] = 1. / (4. / (h ** 2) * np.sin((np.pi * (k + 1) * h / 2.))**2  \
                                                     + 4. / (h ** 2) * np.sin((np.pi * (m + 1) * h / 2.))**2)
 
v = idst(f_h, type=1) * np.sqrt(1 / (2 * N))
p = np.dot(D, v)
solution = dst(p, type=1) * np.sqrt(1 / (2 * N))

true_sol = np.zeros(N)

for k in range(n - 1):
    for m in range(n - 1):
        true_sol[k * (n - 1) + m] = u((k + 1) * h, (m + 1) * h)

print("solution found: ", solution)
print("true solution: ", true_sol)
print("difference: ", np.linalg.norm(solution - true_sol))