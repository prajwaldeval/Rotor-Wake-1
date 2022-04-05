import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = np.array(pd.read_excel('polar DU95W180.xlsx'))[3:,:]
data = np.array(data,dtype='float64')
alpha = data[:,0]
cl = data[:,1]
cd = data[:,2]
clcd = cl/cd
opt = np.max(clcd)
cl_opt = cl[clcd==opt]
cd_opt = cd[clcd==opt]
alpha_opt = alpha[clcd==opt]


plt.plot(alpha,clcd)
plt.grid()
plt.xlabel(r'$\alpha [^\circ]$')
plt.ylabel(r'$C_{l}/C_{d}$')