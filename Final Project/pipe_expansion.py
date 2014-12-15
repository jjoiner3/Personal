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


in2mm = 25.4
OD = 18*in2mm # mm
OD = OD/1000. # m

Pipe_thick = 0.375*in2mm/1000

ID = OD-Pipe_thick*2


nr = 51.

nt = 40000

alpha = 1.172e-5 #thermal diffusity of steel

a = 12.0e-6 #thermal expansion coeff of steel
ft2m = 0.3048
L0 = 200.0 #ft
L0 = L0*ft2m  #m



L = numpy.ones(nt)*L0


dr = Pipe_thick/(nr-1)

Tins = 10.
Ti = numpy.ones(nr)*(Tins-32)*5/9
Topp = 350. #deg F
Ti[0] = (Topp-32)*5/9


sigma = 1/2.0
dt = sigma * dr*dr/alpha



T = Ti.copy()
Lexp2 = L.copy()
for n in range(nt):  

    
    Tn = T.copy() 
    T[1:-1] = Tn[1:-1] + alpha*dt/dr**2*(Tn[2:] -2*Tn[1:-1] + Tn[0:-2])
    T[-1] = T[-2]
    n_stop = n + 1
    
    T_Break = (Topp-32)*5/9
    T_current = T[nr-1]
    
    T_fin = T_Break-T_current
    
    T1 = numpy.mean(Tn)
    T2 = numpy.mean(T)
    dT = T2-T1
    dl = a*L0*dT
    
    Lexp2[n+1] = Lexp[n] + dl
    
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

Lexp = Lexp2/ft2m
print (Lexp2[n_stop]-L[n_stop]/ft2m)*12

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
plt.plot(numpy.linspace(0,time,n_stop), Lexp2[:n_stop], color='#003366', ls='-', lw=3)
plt.ylim(Lexp2[0]*.999,Lex2p[n_stop]*1.001)

locs,labels = plt.yticks()
plt.yticks(locs, map(lambda x: "%.3f" % x, locs))

plt.title('Pipe Length In Time\n')
plt.ylabel('Pipe Length, ft')
plt.xlabel('Time, sec');

plt.tight_layout()
plt.show()

