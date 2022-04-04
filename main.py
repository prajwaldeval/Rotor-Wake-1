import matplotlib.pyplot as plt
import numpy as np
from defs import *

r_cos,dr_cos = geometry_cosine()
r_con,dr_con = geometry_constant()

c = chord(r_con/R)
b = twist(r_con/R)


for i in range(len(r_con)):
    a,ap,Fax,Faz,gamma = solve_anul(U0,r_con[i],dr_con[i],c[i],b[i])
    print(i,'finished with a =',a)
