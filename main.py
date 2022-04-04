import matplotlib.pyplot as plt
import numpy as np
from defs import *


constant = True # True is we want constant, False if we want cosine

r_cos,dr_cos = geometry_cosine(n)
r_con,dr_con = geometry_constant(n)

if constant:
    c = chord(r_con / R)
    b = twist(r_con / R)
    r_arr = r_con
    dr_arr = dr_con
else:
    c = chord(r_cos / R)
    b = twist(r_cos / R)
    r_arr = r_cos
    dr_arr = dr_cos

alphas = np.zeros(n)
phis = np.zeros(n)
a = np.zeros(n)
ap = np.zeros(n)
Fax = np.zeros(n)
Faz = np.zeros(n)
gamma = np.zeros(n)

CT = 0
CP = 0

for i in range(len(r_arr)):
    a[i],ap[i],Fax[i],Faz[i],gamma[i], phi, alpha = solve_anul(U0,r_arr[i],dr_arr[i],c[i],b[i],Prandtl=True, Glauert=True)
    alphas[i] = alpha
    phis[i] = phi
    CT += dr_arr[i] * Fax[i] * nb/(0.5*rho*U0**2*np.pi*R**2)
    CP += dr_arr[i] * Faz[i] * nb * R * omega/(0.5*rho*U0**3*np.pi*R**2)
    print(i,'finished')

print('Finished with CT =',CT,'and CP =',CP)

r_R = r_arr/R

plt.figure()
plt.plot(r_R[1:],a[1:],'r',label='Axial induction factor')
plt.plot(r_R[1:],ap[1:],'b',label='Azimuthal induction factor')

plt.title('TSR = 10')
plt.xlabel('Blade location r/R')
plt.ylabel('Factor [-]')
plt.legend(loc='upper right')
plt.grid()

plt.show()

