# A Library for handling various Methods Of Integration, Sudocode from Fundamentals of System Biology from Marcus Covert.
import numpy as np 

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

#Alec's code

def eulerApprox(y0, t0, tF, step = 1000):
  def diffEq(x, y, deltaT):
    '''
    This is where you want to put the expression you want to find an Euler approximation for (e.g. dy/dx = y)
    '''
    return y + y*deltaT #For dy/dx = y. You can change this
  life = tF - t0
  deltaT = life/step
  ts = np.linspace(t0, tF, step)
  ys = []
  ys.append(y0) #Add initial condition
  for i in range(0, step - 1):
    ys.append(diffEq(ts[i], ys[i], deltaT))
  ys = np.array(ys)
  return (ts, ys)


def runKutApprox(y0, t0, tF, step=1000):
  def diffEq(x, y, deltaT=1):
    '''
    This is where you want to put the expression you want to find an Euler approximation for (e.g. dy/dx = y)
    '''
    return y + y*deltaT #For dy/dx = y. You can change this
  life = tF - t0
  deltaT = life/step
  ts = np.linspace(t0, tF, step)
  ys = []
  ys.append(y0) #Add initial condition
  for i in range(0, step - 1):
      k_1 = diffEq(ts[i], ys[i], 1)
      k_2 = diffEq(ts[i] + (deltaT/2), ys[i] + deltaT*(k_1/2))
      k_3 = diffEq(ts[i] + (deltaT/2), ys[i] + deltaT*(k_2/2))
      k_4 = diffEq(ts[i] + deltaT, ys[i] + deltaT*k_3)

      ys.append(ys[i]+ 1/6 * (k_1 + 2*k_2 + 2*k_3 + k_4) * deltaT) 
  ys = np.array(ys)
  return (ts, ys)


# a test - Euler Method
mode = 2 #Runge-Kutta Test
if mode == 1:
    x, y = eulerApprox(1, 0, 1, step=1000)
    plt.plot(x, y)
    x, y = eulerApprox(1, 0, 1, step=100)
    plt.plot(x, y)
    x, y = eulerApprox(1, 0, 1, step=10)
if mode == 2:
    x, y = runKutApprox(1, 0, 1, step=1000)
    plt.plot(x, y)
    x, y = runKutApprox(1, 0, 1, step=100)
    plt.plot(x, y)
    x, y = runKutApprox(1, 0, 1, step=10)
    plt.plot(x, y)
if mode == 3:
    x, y = alecEuler(1, 0, 1, step=10)
    plt.plot(x, y, label="Euler")
    x, y = runKutApprox(1, 0, 1, step=10)
    plt.plot(x, y, label="Runge-Kutta")
    plt.legend()
plt.show()

while True:
    a = 1
