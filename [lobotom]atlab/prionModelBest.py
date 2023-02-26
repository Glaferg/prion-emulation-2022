# A Library for handling various Methods Of Integration, Sudocode from Fundamentals of System Biology from Marcus Covert.
import numpy as np 

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

#Alec's code

class oilerConfig:
    #Healthy PROT - B
    bDNA = 1
    k_trb = 15
    k_txb = 15
    d_bRNA = 0.002
    d_bPROT = 0.001

    #Prion - B'
    bpDNA = 1
    k_trbp = 15
    k_txbp = 15
    d_bpRNA = 0.002
    d_bpPROT = 0.001

    #Rate of B + B' --> 2B
    k_r = 0.3

def oilerMassKineticApprox(bRNA_i, bPROT_i, bpRNA_i, bpPROT_i, t0, tF, step=1000):
  def diffEq(brna, bprot, bprna, bpprot): #Undercase = mutable values, read massKineticEquations.py to understand this.
    bDNA = oilerConfig.bDNA
    k_trb = oilerConfig.k_trb
    k_txb = oilerConfig.k_txb
    d_bRNA = oilerConfig.d_bRNA
    d_bPROT = oilerConfig.d_bPROT

    bpDNA = oilerConfig.bpDNA
    k_trbp = oilerConfig.k_bptr
    k_txbp = oilerConfig.k_txbp
    d_bpRNA = oilerConfig.d_bpRNA
    d_bpPROT = oilerConfig.d_bpPROT

    '''
    This is where you want to put the expression you want to find an Euler approximation for (e.g. dy/dx = y)
    '''

    bRNA = k_txb*bDNA - (d_bRNA * brna)
    bPROT = k_trb*brna - k_r*bprot*bpprot - (d_bPROT * bprot))

    bpRNA = k_txbp*bpDNA - (d_bpRNA * bprna)
    bpPROT = k_trbp*bprna + k_r*bprot*bpprot - (d_bpPROT * bpprot))

    return (bRNA, bPROT, bpRNA, bpPROT) #For dy/dx = y. You can change this
  life = tF - t0
  deltaT = life/step
  ts = np.linspace(t0, tF, step)
  bRNAs = []
  bpRNAs = []
  bPROTs = []
  bpPROTs = []
  bRNAs.append(bRNA_i) #Add initial condition
  bpRNAs.append(bpRNA_i)
  bPROTs.append(bPROT_i)
  bpPROTs.append(bpPROT_i)
  for i in range(0, step - 1):
      k_1 = diffEq(bRNAs[i], bPROTs[i], bpRNAs[i], bpPROTs[i])
      # k_2 = diffEq(RNAs[i] + (deltaT/2), ys[i] + deltaT*(k_1/2))
      # k_3 = diffEq(ts[i] + (deltaT/2), ys[i] + deltaT*(k_2/2))
      # k_4 = diffEq(ts[i] + deltaT, ys[i] + deltaT*k_3)
      bRNAs.append(k_1[0]*deltaT + bRNAs[i])
      bPROTs.append(k_1[1]*deltaT + bPROTs[i])
      bpRNAs.append(k_1[2]*deltaT + bpRNAs[i])
      bpPROTs.append(k_1[3]*deltaT + bpPROTs[i])

  bRNAs = np.array(bRNAs)
  bPROTs = np.array(bPROTs)
  bpRNAs = np.array(bpRNAs)
  bpPROTs = np.array(bpPROTs)

  return [bRNAs, bPROTs, ts, bpRNAs, bpPROTs]

rna_i = 100000
prot_i = 1
results = oilerMassKineticApprox(bRNA_i=rna_i, bPROT_i=prot_i, bpRNA_i=rna_i, bpPROT_i=prot_i, t0=0, tF=10000, step=1000090)
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
