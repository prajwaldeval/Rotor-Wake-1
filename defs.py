import numpy as np
import matplotlib.pyplot as plt

def chord(mu):
    return 0.18 - 0.06*mu

def twist(mu):
    return -50*mu + 35

def geometry_cosine():
    r = np.linspace(0,np.pi,n)
    r = (1-np.cos(r))/2
    blen = R-r_hub
    r = r_hub + blen*r
    return r

def geometry_constant():
    r = np.linspace(r_hub, R, n)
    return  r



def aero_coeffs(alpha,filename = 'ARAD8pct_polar.txt'):
    data = np.genfromtxt(filename,skip_header=2)
    cl = np.interp(alpha,data[:,0].data[:,1])
    cd = np.interp(alpha,data[:,0].data[:,2])
    return cl, cd


def BE_loads(a,ap,r,b,c):
    Vtan = RPM*r*(1+ap)
    Vax = U0*(1-a)
    Vps= Vtan**2+Vax**2
    phi = np.atan(Vax,Vtan)
    alpha = phi - b
    cl, cd = aero_coeffs(alpha)
    L = 0.5*c*rho*Vps*cl
    D = 0.5*c*rho*Vps*cd
    Faz = L*np.sin(phi) - D*np.cos(phi)
    Fax = L*np.cos(phi) + D*np.sin(phi)

    return Fax,Faz

def MT_induction(Fax,Faz,r,b,c,Glauert):
    # CT = ()/()
    return a,ap




if __name__=='__main__':
        R = 0.7 #m
        r_hub = 0.25*R
        b_col = 46 #deg
        U0 = 60 #m/s
        RPM = 1200
        h = 2000 #m
        n = 50 #number of discretised blade elements
        r_arr_cos = geometry_cosine()
        r_arr_con = geometry_constant()
        aero = np.genfromtxt('ARAD8pct_polar.txt',skip_header=2)
        plt.scatter(aero[:,0],aero[:,1])

