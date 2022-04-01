import matplotlib.pyplot as plt
import numpy as np
from defs import *

r_cos,dr_cos = geometry_cosine()
r_con,dr_con = geometry_constant()

c = chord(r_cos/R)
b = twist(r_cos/R)


for i in range(len(r_cos)):
    a,ap,Fax,Faz,gamma = solve_anul(U0,r_cos[i],dr_cos[i],c[i],b[i])
    print(a,ap)