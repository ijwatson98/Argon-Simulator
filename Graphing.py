# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 21:42:37 2021

@author: Surface
"""

import matplotlib.pyplot as plt


X1, Y1 = [], []
for line in open('melting_obs1.2.txt', 'r'):
  values = [float(s) for s in line.split()]
  X1.append(values[0])
  Y1.append(values[4])
  
X2, Y2 = [], []
for line in open('melting_obs1.21.txt', 'r'):
  values = [float(s) for s in line.split()]
  X2.append(values[0])
  Y2.append(values[4])
  
X3, Y3 = [], []
for line in open('melting_obs1.22.txt', 'r'):
  values = [float(s) for s in line.split()]
  X3.append(values[0])
  Y3.append(values[4])


plt.title('Argon: MSD Comparison')
plt.xlabel('time/tau')
plt.ylabel('r/sigma^2')
plt.plot(X1, Y1, color = 'green')
plt.plot(X2, Y2, color = 'blue')
plt.plot(X3, Y3, color = 'red')
plt.legend(['T = 1.20', 'T = 1.21', 'T = 1.22'])
plt.show()
