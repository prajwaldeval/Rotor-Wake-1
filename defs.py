import numpy as np
import matplotlib.pyplot as plt

R = 0.7 #m
r_hub = 0.25*R
b_col = 46 #deg
U0 = 60 #m/s
RPM = 1200
omega=RPM*2*np.pi/60 #rad/s
h = 2000 #m
nb = 3 #number of blades
n = 50 #number of discretised blade elements
rho = 1.00649 #kg/m3 at 2000m altitude (ISA)

def chord(mu):
    return 0.18 - 0.06*mu

def twist(mu):
    return -50*mu + 35

def geometry_cosine():
    r = np.linspace(0,np.pi,n+1)
    r = (1-np.cos(r))/2
    blen = R-r_hub
    r = r_hub + blen*r
    dr = r[1:]-r[:-1]
    r = (r[1:]+r[:-1])/2
    return r, dr

def geometry_constant():
    r = np.linspace(r_hub, R, n+1)
    dr = r[1:] - r[:-1]
    r = (r[1:] + r[:-1]) / 2
    return  r,dr



def aero_coeffs(alpha,filename = 'ARAD8pct_polar.txt'):
    data = np.genfromtxt(filename,skip_header=2)
    cl = np.interp(alpha,data[:,0].data[:,1])
    cd = np.interp(alpha,data[:,0].data[:,2])
    return cl, cd

def BE_loads(a,ap,r,dr,b,c):
    Vtan = RPM*r*(1+ap)
    Vax = U0*(1+a) #1-a
    Vps= Vtan**2+Vax**2
    phi = np.atan(Vax,Vtan)
    b = b + b_col
    alpha = - phi + b
    cl, cd = aero_coeffs(alpha)
    L = 0.5*c*rho*Vps*cl
    D = 0.5*c*rho*Vps*cd
    Faz = L*np.sin(phi) - D*np.cos(phi)
    Fax = L*np.cos(phi) + D*np.sin(phi)
    gamma = 0.5 * np.sqrt(Vps) * cl * c  # vorticity
    return Vax,Vtan,Fax,Faz,gamma

def MT_induction(Fax,Faz,r,dr,b,c,Glauert,Prandtl):
    CT = (Fax*nb*r)/(0.5*rho*U0**2*2*np.pi*r*dr)
    ap = (Faz*nb)/(2*rho*(2*np.pi*r)*U0*(1+a)*r*omega) #1-a

    a = 0.5 - 0.5*np.sqrt(1-CT)
    if Glauert:
        CT1 = 1.816
        CT2 = 2*np.sqrt(CT1) - CT1
        if CT>=CT2:
            a = 1+ (CT-CT1)/(4*np.sqrt(CT1)-4)

    if Prandlt:
        mu = r/R
        l = omega*R/U0
        p1 = -nb/2*(1-mu)/(mu)*np.sqrt(1+(l**2*mu**2)/(1-a)**2)
        ftip = 2/np.pi * np.arccos(np.exp(p1))
        p2 = -nb/2*(mu-r_hub/R)/(mu)*np.sqrt(1+(l**2*mu**2)/(1-a)**2)
        froot = 2/np.pi * np.arccos(np.exp(p2))
        ftot = ftip + froot
        a = a/ftot
        ap = ap/ftot

    return a,ap

def solve_anul(Uinf, r, dr, c, b):

    #initialization of variables
    a = 0   # axial induction factor
    ap = 0  # tangential induction factor

    convergence_error = 0.00001
    max_nr_iterations = 100

    for i in range(max_nr_iterations):

        Vax, Vtan, Fax, Faz, gamma = BE_loads(a,ap,r,dr,b,c)
        an, apn = MT_induction(Fax,Faz,r,dr,b,c, Prandtl=True, Glauert=True)

        anext = 0.25*a + 0.75*an
        apnext = 0.25*ap + 0.75*apn

        if (np.abs(a-anext)<convergence_error):
            print(i)
            break
        if i==max_nr_iterations:
            print('not conveerged')
        a = anext
        ap = apnext

    return a,ap,Fax,Faz,gamma


if __name__=='__main__':
        R = 0.7 #m
        r_hub = 0.25*R
        b_col = 46 #deg
        U0 = 60 #m/s
        RPM = 1200
        omega=RPM*2*np.pi/60 #rad/s
        h = 2000 #m
        nb = 3 #number of blades
        n = 50 #number of discretised blade elements
        rho = 1.00649 #kg/m3 at 2000m altitude (ISA)
        r_arr_cos,dr_cos = geometry_cosine()
        r_arr_con,dr_con = geometry_constant()


