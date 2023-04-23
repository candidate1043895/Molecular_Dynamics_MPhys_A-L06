import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.stats import boltzmann
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, k
from math import exp as exp
from mpmath import besselk
import mpmath as mp
from scipy.constants import physical_constants
from scipy.special import kn
import math
from math import pi
from scipy.optimize import curve_fit
from scipy.constants import epsilon_0
mp.pretty=True
u=physical_constants['atomic mass constant'][0]
e=physical_constants['atomic unit of charge'][0]
kB=physical_constants['Boltzmann constant in eV/K'][0]
maxwell = stats.maxwell



#file location, coord.txt contains
#           'some nonsense strings we want to get rid of'
#           ATOMS_id type x y z vx vy vz c_px c_py c_pz
#           data
where_your_coord_file_lives='/Users/candidate123456/MPhys_MD/lammps-rel/examples/USER/rel/magneto-statics/18_runs_e/ewald'
coord_file='{}/coord.txt'
file=coord_file.format(where_your_coord_file_lives)

#find number of particles in the system
def particle_number(File):
    f=open(File)
    lines=f.readlines()
    line=lines[3]
    n=int(line.split()[0])
    return n
 
 
def velocity_range0(File, lower, dump, pn):
    
    testFile = File
    #file = open(testFile, 'r')
    lowindex=(lower/dump)*(pn+9)
    data=[]
    with open(testFile, 'r') as f:
        lines=f.readlines()[int(lowindex):]
             
        for row in lines:
            if 'ITEM' not in row:
                new_row=re.split(r'\s', row)
                if len(new_row)==13:
                    new_row2=list(filter(None, new_row))
                    new_row3=[float(n) for n in new_row2]
                    data.append(new_row3)
    
    data=np.array(data)
    vx=data[:,5]
    vy=data[:,6]
    vz=data[:,7]
    v=[]
    for i in range(0, len(data), 1):
        vtot=100*(vx[i]**2+vy[i]**2+vz[i]**2)**0.5
        v.append(vtot/(3e8))
        
    return v




def etot_range(File, lower, pn):
    
    testFile = File
    lowindex=(lower/10000)*(pn+9)
    data=[]
    with open(testFile, 'r') as f:
        lines=f.readlines()[int(lowindex):]
             
        for row in lines:
            if 'ITEM' not in row:
                new_row=re.split(r'\s', row)
                if len(new_row)==13:
                    new_row2=list(filter(None, new_row))
                    new_row3=[float(n) for n in new_row2]
                    data.append(new_row3)
    
    data=np.array(data)
    vx=data[:,5]
    vy=data[:,6]
    vz=data[:,7]
    v=[]
    for i in range(0, len(data), 1):
        vtot=100*(vx[i]**2+vy[i]**2+vz[i]**2)**0.5
        v.append(vtot/(3e8))
        
    return v

def velocity(File):
    testFile = File
    file = open(testFile, 'r')
    data=[]
    for row in file:
        if 'ITEM' not in row:
            new_row=re.split(r'\s', row)
            if len(new_row)==13:
                new_row2=list(filter(None, new_row))
                new_row3=[float(n) for n in new_row2]
                data.append(new_row3)
    
    data=np.array(data)
    vx=data[:,5]
    vy=data[:,6]
    vz=data[:,7]
    v=[]
    for i in range(0, len(data), 1):
        vtot=100*(vx[i]**2+vy[i]**2+vz[i]**2)**0.5
        v.append(vtot/(3e8))
        
    return v

def cenbin(bins):
    cen_bin=[]
    for i in range(0, len(bins)-1, 1):
        binv=(bins[i]+bins[i+1])/2
        cen_bin.append(binv)
    return cen_bin

def plot_hist(num_bins, x_max, v, yscale, label, name):
    n, bins, patches = plt.hist(v, num_bins, density = True,
                                alpha = 0.7, histtype='step', label=label)
    plt.xlabel('v/c')
    plt.ylabel('f(v)')
    plt.xlim([0, x_max])
    if yscale=='log':
        plt.yscale('log', nonposy='clip')
        
    plt.title(name,fontweight ="bold")
    cen_bin=cenbin(bins)
    
    return n, cen_bin

def moving_avg(x, y, n):
    n=int(n)
    xm=[]
    ym=[]
    for i in range(0, len(x)+1-n, 1):
        x0=sum(x[i:i+n])/n
        y0=sum(y[i:i+n])/n
        xm.append(x0)
        ym.append(y0)
    return xm, ym
    
    
def boltz_fit(v,xmax, label, lw):
    params = maxwell.fit(v, floc=0)

    x = np.linspace(0, xmax, 100)
    #plt.plot(cenb_12, n_12, label=label)
    plt.plot(x, maxwell.pdf(x, *params), lw=lw, label='Bfit_{}'.format(label))
    plt.legend(loc='upper right')
    print(label, params)
    

def MJ_dist(mass, temp, n, beta):
    m=mass
    T=temp
    theta=(m*c**2)/(k*T)
    K_n=float(kn(n, theta))
    gamma=1/(1-(beta)**2)**0.5
    MJ=((beta**2*gamma**5)/(K_n/theta))*(np.exp(-theta*gamma))
    
    #source:  https://en.wikipedia.org/wiki/Maxwell–Jüttner_distribution
    
    return MJ, theta, K_n, gamma

#MJ distribution for Ar
def MJ_function(beta, temp):
    m=40*u
    T=temp
    theta=(m*c**2)/(k*T)
    K_n=float(kn(2, theta))
    gamma=1/(1-(beta)**2)**0.5
    MJ=((beta**2*gamma**5)/(K_n/theta))*(np.exp(-theta*gamma))
    
    #source:  https://en.wikipedia.org/wiki/Maxwell–Jüttner_distribution
    
    return MJ

#MJ distribution for electrons
def MJ_function_e(beta, temp):
    m=0.00055*u
    T=temp
    theta=(m*c**2)/(k*T)
    K_n=float(kn(2, theta))
    gamma=1/(1-(beta)**2)**0.5
    MJ=((beta**2*gamma**5)/(K_n/theta))*(np.exp(-theta*gamma))
    
    #source:  https://en.wikipedia.org/wiki/Maxwell–Jüttner_distribution
    
    return MJ




pnm=particle_number(file)
v=velocity_range0(file, 0.2e7, 100000, pnm)
nlj_bins=100
x_max=0.15
n, cb=plot_hist(nlj_bins, x_max, v, 'linear', 'N/A', 'ELECTRON coul')

p0=[6e6]
popt, pcov = curve_fit(MJ_function_e, cb, n, p0 )
print('fitted temperature',popt, pcov)
print('percentage error', np.sqrt(pcov.diagonal())/popt[0])
plt.figure()
plt.plot(cb, n, '.', label='exp data')
beta=np.linspace(0, x_max, 500)
plt.plot(beta, MJ_function_e(np.array(beta),*popt),color='red', label='fit_T={}'.format(f"{popt[0]:.7}"))
#params = maxwell.fit(vew_new2_mag, floc=0)
#plt.plot(beta, maxwell.pdf(beta, *params), lw=2, color='red', label='Boltzmann dist')
print('maxwell parameters',params )
#plt.plot(beta, y)
#plt.yscale('log')
plt.title('what_ever_you_want_to_name_your_plot',fontweight ="bold")
plt.legend()
plt.xlim(0)
plt.plot()
plt.ylabel('f(v)')
plt.xlabel('v/c')
plt.savefig('what_ever_you_want_to_name_your_plot.jpg', format = 'jpg',bbox_inches='tight', dpi=600)
