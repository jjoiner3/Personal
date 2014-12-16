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
alpha = 1.172e-5 #thermal diffusity, m²/s of steel
a = 12.0e-6 #thermal expansion coeff, 10-6/K, of steel
nu = 0.3 # Poisson's ratio streel, carbon steel
E = 205.5e9*0.000145 # Modulous of Elasticity, Pa, carbon steel

c13 = E*(nu)/((1.+nu)*(1.-2.*nu))
c12 =c13
c23 = c13

c33 = E*(1.-nu)/((1.+nu)*(1.-2.*nu))
Beta = (c13+c23+c33)*a
#Beta = (c33)*a
2022169536
A = c33*a+Beta
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

#dt= sigma * dr*dr/A/alpha
#dt2= sigma * dr*dr/(alpha**2)
#dt3= sigma * dr*dr/alpha


dt = sigma * dr*dr/alpha


# Initial Conditions
Tins = 10. # Install Temperature deg F
Tins = ((Tins-32)*5/9 + 274.15) #converted to Kelvin
Ti = numpy.ones((nr,nt))*Tins#*Tins#*((Tins-32)*5/9 + 274.15)# Initial Temperature array converted to deg C
Topp = 350. # Operational Temperature deg F
Ti[0,:] = (Topp-32)*5/9  + 274.15#Sets the temperature at the inner diameter to Operational Temp, deg C  
L = numpy.ones((nr,nt))*L0 #Initalizes array for lengths
#dL = numpy.zeros(nr)

# Final Temperature and Pipe Lengths
T = Ti.copy()
Lexp = L.copy()

# Stress and Strain Arrays, Initailly Zero Stress or Strain
StressI = numpy.zeros((nr,nt))
StrainI = numpy.zeros((nr,nt))

# Final Stress and Strain
StressF = StressI.copy()
StrainF = StressI.copy()

# Final Temperature
T_Break = (Topp-32)*5/9 + 274.15

for n in range(1,nt):  

#   Solve for Temp
    Tn = T.copy() 
    T[1:-1,n] = Tn[1:-1,n-1] + alpha*dt/dr**2*(Tn[2:,n-1] -2*Tn[1:-1,n-1] + Tn[0:-2,n-1])
    T[-1,n] = T[-2,n]
    
    #StressN = StressF.copy()
    #StressF[1:-1] = StressN[1:-1] + A*alpha*dt/dr**2*(Tn[2:] -2*Tn[1:-1] + Tn[0:-2])
    #StressF[0] = StressF[1]
    #StressF[-1] = StressF[-2]


    StrainN = StrainF.copy()
    StrainF[1:-1,n] = StrainN[1:-1,n-1] + (alpha**2)*dt/dr**2*(Tn[2:,n-1] -2*Tn[1:-1,n-1] + Tn[0:-2,n-1])
    StrainF[0,n] = StrainF[1,n]
    StrainF[-1,n] = StrainF[-2,n]

    Lexp[:,n] = StrainF[:,n]*L[:,n]+L[:,n]
    
    #Lexp = Lexp+dL
    
# Determining Stop Condition
    T_current = T[nr-1,n] #Temp at OD
    T_fin = T_Break-T_current #Difference between Temp at OD and Temp at ID
    n_stop = n + 1
    if Topp > Tins:
        if T_fin <= 0.1:
            print ("Outside of pipe reached 350F at time {0:.2f}s, in time step {1:d}.".format(dt*n_stop, n_stop))
            break
    else: 
        if T_fin >= -0.1:
            print ("Outside of pipe reached 350F at time {0:.2f}s, in time step {1:d}.".format(dt*n_stop, n_stop))
            break         


StressF3 = E*StrainF #psi

T = ((T - 274.15)*9./5.)+32.
time = n*dt

#print n_stop, time, T[50,nt]>=Topp*0.1

Lexp = Lexp/ft2m
#print (Lexp[nr-1,nt]-L[nr-1,nt]/ft2m)*12


