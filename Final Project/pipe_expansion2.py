# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy
import matplotlib.pyplot as plt
#matplotlib inline
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16




#Material Properties
alpha = 1.172e-5 #thermal diffusity of steel
a = 12.0e-6 #thermal expansion coeff of steel
nu = 0.29 # Poisson's ratio streel, carbon steel
E = 190.0e6 # Modulous of Elasticity GPa, carbon steel

c13 = E*(nu)/((1.+nu)*(1.-2.*nu))
c12 =c13
c23 = c13

c33 = E*(1.-nu)/((1.+nu)*(1.-2.*nu))
Beta = (c13+c23+c33)*alpha


A = c33*a*alpha **2+Beta
# Conversion Factors
in2mm = 25.4 #converts to SI units
ft2m = 0.3048 #Converts ft to meters

# Pipe Size
OD = 18 # Outer Diameter in inches
OD = OD*in2mm/1000. # converts inches to meters
Pipe_thick = 0.375*in2mm/1000 # converts inches to meters
ID = OD-Pipe_thick*2 # Inner Diameter

#Pipe Length
L0 = 200.0 # Initial length of pipe, ft
L0 = L0*ft2m  #Converts ft to meters


# Grid Size
nr = 51.
nt = 40000
dr = Pipe_thick/(nr-1)
# Stability Requirements
sigma = 1/2.0
dt = sigma * dr*dr/alpha


# Initial Conditions
Tins = 95. # Install Temperature deg F
Ti = numpy.ones(nr)*(Tins-32)*5/9 # Initial Temperature array converted to deg C
Topp = 42. # Operational Temperature deg F
Ti[0] = (Topp-32)*5/9 #Sets the temperature at the inner diameter to Operational Temp, deg C  
L = numpy.ones(nt)*L0 #Initalizes array for lengths at each time step


# Final Temperature and Pipe Lengths
T = Ti.copy()
Lexp = L.copy()

# Stress and Strain Arrays, Initailly Zero Stress or Strain
StressI = numpy.zeros(nr)
StrainI = numpy.zeros(nr)

# Final Stress and Strain
StressF = StressI.copy()
StrainF = StressI.copy()

# Final Temperature
T_Break = (Topp-32)*5/9

for n in range(nt):  

#   Solve for Temp
    Tn = T.copy() 
    T[1:-1] = Tn[1:-1] + alpha*dt/dr**2*(Tn[2:] -2*Tn[1:-1] + Tn[0:-2])
    T[-1] = T[-2]
    
    StressN = StressF.copy()
    StressF[1:-1] = StressN[1:-1] + A*dt/dr**2*(Tn[2:] -2*Tn[1:-1] + Tn[0:-2])
    StressF[-1] = StressF[-2]


    StrainN = StrainF.copy()
    StrainF[1:-1] = StrainN[1:-1] + (alpha**2)*dt/dr**2*(Tn[2:] -2*Tn[1:-1] + Tn[0:-2])
    StrainF[-1] = StrainF[-2]


# Determining Stop Condition
    T_current = T[nr-1] #Temp at OD
    T_fin = T_Break-T_current #Difference between Temp at OD and Temp at ID
    n_stop = n + 1
    if Topp > Tins:
        if T_fin <= 0.001:
            print ("Outside of pipe reached 42F at time {0:.2f}s, in time step {1:d}.".format(dt*n_stop, n_stop))
            break
    else: 
        if T_fin >= -0.001:
            print ("Outside of pipe reached 42F at time {0:.2f}s, in time step {1:d}.".format(dt*n_stop, n_stop))
            break         
         
    
T = (T*9./5.)+32.
time = n*dt

print n_stop, time, T[50]>=Topp*0.1

Lexp = Lexp/ft2m
print (Lexp[n_stop]-L[n_stop]/ft2m)*12

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

plt.rc('font', **font)
plt.figure(num=None, figsize=(10, 8), dpi=80, facecolor='w', edgecolor='k')


plt.figure(1)
plt.subplot(212)
plt.plot(numpy.linspace(ID*1000/in2mm,OD*1000/in2mm,nr), T, color='#003366', ls='-', lw=3)
plt.ylim(0,Topp*1.2)

plt.title('Pipe Cross Section Temp at Final Time\n')
plt.xlabel('Thickness Of Pipe, in')
plt.ylabel('Temperature, F');


#Lexp_in[:n_stop] = (Lexp[:n_stop] - (Lexp[:n_stop]))*12/ft2m
plt.subplot(211)
plt.plot(numpy.linspace(0,time,n_stop), Lexp[:n_stop], color='#003366', ls='-', lw=3)
plt.ylim(Lexp[0]*.999,Lexp[n_stop]*1.001)

locs,labels = plt.yticks()
plt.yticks(locs, map(lambda x: "%.3f" % x, locs))

plt.title('Pipe Length In Time\n')
plt.ylabel('Pipe Length, ft')
plt.xlabel('Time, sec');

plt.tight_layout()
plt.show()

