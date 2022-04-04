import matplotlib.pyplot as plt
import numpy as np
from defs import *

n=50
r_cos,dr_cos = geometry_cosine(n)
r_con,dr_con = geometry_constant(n)

c = chord(r_cos/R)
b = twist(r_cos/R)

r_arr = r_cos
a = np.zeros(n)
ap = np.zeros(n)
Fax = np.zeros(n)
Faz = np.zeros(n)
gamma = np.zeros(n)

CT = 0
CP = 0

for i in range(len(r_cos)):
    a[i],ap[i],Fax[i],Faz[i],gamma[i] = solve_anul(U0,r_con[i],dr_con[i],c[i],b[i],Prandtl=True, Glauert=True)
    CT += dr_cos[i] * Fax[i] * nb/(0.5*rho*U0**2*np.pi*R**2)
    CP += dr_cos[i] * Faz[i] * nb * R * omega/(0.5*rho*U0**3*np.pi*R**2)
    print(i,'finished')

print('Finished with CT =',CT,'and CP =',CP)
# CT = np.sum(dr*results[:,3]*NBlades/(0.5*rho*U0**2*np.pi*R**2))
# CP = np.sum(dr*results[:,4]*results[:,2]*NBlades*Radius*Omega/(0.5*Uinf**3*np.pi*Radius**2))

plt.plot(r_arr,a,'r')
plt.plot(r_arr,ap,'b')

plt.show()

