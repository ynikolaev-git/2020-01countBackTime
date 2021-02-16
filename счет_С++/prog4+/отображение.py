# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 08:56:07 2020

@author: 1
"""

## ==========================================================
##     рисуем скорость звука
#import pandas as pd
#import matplotlib.pyplot as plt
#net = pd.read_csv('net.csv', sep=';')
#piston = pd.read_csv('piston.csv', sep=';')
#
#rw = 1.07353
#pNum = 182
#rmin = float(piston.query('num == ' + str(pNum)).values[:,2])
#print('rmin = ', rmin)
#t = float(piston.query('num == ' + str(pNum)).values[:,1])
#rmax = rw
#print('rmax = ', rmax)
#nr = 100
#dr = (rmax-rmin)/nr
#i = 1
#r = rmin + i*dr
#
#for i in range(nr):
#    r = rmin + i*dr
#    net['met'] = (net['t']-t)**2 + (net['x']-r)**2
#    net = net.sort_values(by=['met'])
#    c = net.iat[0,5]
#    plt.plot(r, c, 'ro')
#plt.axis([rmin-rmin/10000, 1.0736, 1.4, 1.5])
#plt.show()

# ==========================================================
#  рисуем всю область
import pandas as pd
import matplotlib.pyplot as plt
net = pd.read_csv('regnet2.csv', sep=';')

t = net.values[:, 4]
x = net.values[:, 6]
t_min = min(t)
t_max = max(t)
print('t_min =', t_min)
print('t_max =', t_max)
x_min = min(x)
x_max = max(x)
print('x_min =', x_min)
print('x_max =', x_max)

net = net.query('numOfLayer == 48')
t = net.values[:, 4]
x = net.values[:, 6]

#x_min = min(x)
#print('x_min =', x_min)
#x_max = max(x)
#print('x_max =', x_max)
#plt.figure(figsize=(10, 5))
#plt.axis(xlim=(x_min, x_max), ylim=(t_min, t_max))
plt.plot(t, x, 'ro')
plt.show()

# ==========================================================


#i = 0
#while (i<10):
#    df = table.query('numOfLayer == ' + str(i))
#    i += 1
#    t = df.values[:, 2]
#    x = df.values[:, 3]
##    plt.figure(figsize=(10, 5))
##    plt.plot(x, t, 'ro')
#    plt.plot(x, t)
#    plt.draw()
#    plt.pause(0.001)
##    plt.clf()
#plt.ioff()
#plt.show()


# ==========================================================
#df = table.query('numOfLayer == 350')
#df = net
#df = df.append(table.query('numOfLayer==426'))
#print(df)

#t = df.values[:, 2]
#x = df.values[:, 3]
#plt.plot(x, t, 'ro')
#plt.show()

# ==========================================================


