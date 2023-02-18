# A Library for handling various Methods Of Integration, Sudocode from Fundamentals of System Biology from Marcus Covert.
import numpy as np 

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

#Alec's code

class oilerConfig:
    DNA = 1
    k_tr = 15
    k_tx = 15
    d_RNA = 0.002
    d_PROT = 0.001
def oilerMassKineticApprox(RNA_i, PROT_i, t0, tF, step=1000):
  def diffEq(rna, prot):
    DNA = oilerConfig.DNA
    k_tr = oilerConfig.k_tr
    k_tx = oilerConfig.k_tx
    d_RNA = oilerConfig.d_RNA
    d_PROT = oilerConfig.d_PROT

    '''
    This is where you want to put the expression you want to find an Euler approximation for (e.g. dy/dx = y)
    '''
    RNA = k_tx*DNA - (d_RNA * rna)
    PROT = k_tr*rna - (d_PROT * prot)
    return (RNA, PROT) #For dy/dx = y. You can change this
  life = tF - t0
  deltaT = life/step
  ts = np.linspace(t0, tF, step)
  RNAs = []
  PROTs = []
  RNAs.append(RNA_i) #Add initial condition
  PROTs.append(PROT_i)
  for i in range(0, step - 1):
      k_1 = diffEq(RNAs[i], PROTs[i])
      # k_2 = diffEq(RNAs[i] + (deltaT/2), ys[i] + deltaT*(k_1/2))
      # k_3 = diffEq(ts[i] + (deltaT/2), ys[i] + deltaT*(k_2/2))
      # k_4 = diffEq(ts[i] + deltaT, ys[i] + deltaT*k_3)
      RNAs.append(k_1[0]*deltaT + RNAs[i])
      PROTs.append(k_1[1]*deltaT + PROTs[i])
  RNAs = np.array(RNAs)
  PROTs = np.array(PROTs)
  return [RNAs, PROTs, ts]

rna_i = 100000
prot_i = 1
results = oilerMassKineticApprox(RNA_i=rna_i, PROT_i=prot_i, t0=0, tF=10000, step=1000090)
mode = 3
if mode == 1:
  plt.plot(results[2], results[0]) #RNA/TIME
  plt.title(f"MAC-PROTEIN: RNA/TIME @ {rna_i} RNA, {prot_i} PROT")
  print(f"k_tr = {oilerConfig.k_tr}, k_tx = {oilerConfig.k_tx}, d_rna = {oilerConfig.d_RNA}, d_prot = {oilerConfig.d_PROT}, DNA = {oilerConfig.DNA}")
elif mode == 2:
  plt.plot(results[2], results[1]) #PROT/TIME
  plt.title(f"MAC-PROTEIN: PROT/TIME @ {rna_i} RNA, {prot_i} PROT")
  print(f"k_tr = {oilerConfig.k_tr}, k_tx = {oilerConfig.k_tx}, d_rna = {oilerConfig.d_RNA}, d_prot = {oilerConfig.d_PROT}, DNA = {oilerConfig.DNA}")
else:
  plt.plot(results[2], results[0], label="RNA") #RNA/TIME
  plt.title(f"MAC-PROTEIN: RNA/TIME @ {rna_i} RNA, {prot_i} PROT")
  print(f"k_tr = {oilerConfig.k_tr}, k_tx = {oilerConfig.k_tx}, d_rna = {oilerConfig.d_RNA}, d_prot = {oilerConfig.d_PROT}, DNA = {oilerConfig.DNA}")
  plt.plot(results[2], results[1], label="PROT") #PROT/TIME
  plt.title(f"MAC-PROTEIN: PROT/TIME @ {rna_i} RNA, {prot_i} PROT")
  print(results[1])
  plt.legend()

plt.show()

while True:
    a = 1
