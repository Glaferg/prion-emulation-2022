# A Library for handling various Methods Of Integration, Sudocode from Fundamentals of System Biology from Marcus Covert.

import numpy as np 

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

#Define two arrays to hold the time and concentration data: tarr = [t0]; xarr = [x0];

tarr = np.array(np.empty(2), dtype=np.int64)
xarr = np.array(np.empty(2), dtype=np.int64)

def oilerApprox(tF, xF, function, t0=1, x0=1,  step=1):
    def f(x):
            #use the f-string expression {x} in function-string to substitute in current x
            fun = f"{function}"
            return(eval(f"{fun}"))

    life = tF - t0
    deltaT = life/step

    tarr = np.array([t0], dtype=np.int64)
    xarr = np.array([x0], dtype=np.int64)
    tnew = np.array([t0], dtype=np.int64)
    xnew = np.array([x0], dtype=np.int64)

    for i in range(1, tF):
        tnew = tarr[i - 1] + deltaT
        xnew = xarr[i - 1] + f(xarr[i - 1]) * deltaT
        tnew, xnew = np.array()
        tarr = np.concatenate((tarr, tnew))
        xarr = np.concatenate((xarr, xnew))
    return(tarr, xarr)

def plot(x, y):
    fig, ax = plt.subplots()  # Create a figure containing a single axes.
    ax.plot((x, y));  # Plot some data on the axes.
    plt.show()

# a test
oilerApprox(tF=100, xF=100, function="{x} + 2")
plot(tarr, xarr)
print(tarr)
print(xarr)
while True:
    a = 1
