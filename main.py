import matplotlib.pyplot as plt
import numpy as np
from defs import *

n=50
r_cos,dr_cos = geometry_cosine(n)
r_con,dr_con = geometry_constant(n)

c = chord(r_cos/R)
b = twist(r_cos/R)

r_arr = r_cos
a =  np.zeros(n)
ap = np.zeros(n)
Fax = np.zeros(n)
Faz = np.zeros(n)
gamma = np.zeros(n)

for i in range(len(r_cos)):
    a[i],ap[i],Fax[i],Faz[i],gamma[i] = solve_anul(U0,r_con[i],dr_con[i],c[i],b[i],Prandtl=False, Glauert=True)
    print(i,'finished')
plt.plot(r_arr,ap)
