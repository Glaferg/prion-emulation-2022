# A Model of Prions derived from Alec Lourenco's teaching and the principles of Mass Action Kinetics. Written 2023 by then-Bridges student and Agoura Hills
# Aperture Laboratories Consumer Advocate  Gabriel LaFrance Fergesen.

import numpy as np 

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

#Alec's code

class oilerConfig:
    #Rate of B + B' --> 2B'
    k_r = 0.3

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

def oilerMassKineticApprox(RNA_i, PROT_i, t0, tF, step=1000):
  def diffEq(rna, prot):
    bpDNA = oilerConfig.bpDNA
    k_trbp = oilerConfig.k_bptr
    k_txbp = oilerConfig.k_txbp
    d_bpRNA = oilerConfig.d_bpRNA
    d_bpPROT = oilerConfig.d_bpPROT

    bDNA = oilerConfig.bDNA
    k_trb = oilerConfig.k_trb
    k_txb = oilerConfig.k_txb
    d_bRNA = oilerConfig.d_bRNA
    d_bPROT = oilerConfig.d_bPROT

    '''
    This is where you want to put the expression you want to find an Euler approximation for (e.g. dy/dx = y)
    '''
    bRNA = k_txb*bDNA - (d_bRNA * bRNA)
    bPROT = k_trb*brna - k_r*bPROT*bpPROT - (d_bPROT * bprot)
    
    bpRNA = k_txbp*bpDNA - (d_bpRNA * bpRNA)
    bpPROT = k_trbp*bpRNA + k_r*bPROT*bpPROT - (d_bPROT * bPROT)
    
    return (bRNA, bPROT) #For dy/dx = y. You can change this
  life = tF - t0
  deltaT = life/step
  ts = np.linspace(t0, tF, step) #Generates list of times.
  bRNAs = []
  bPROTs = []
  RNAs.append(RNA_i) #Add initial condition
  PROTs.append(PROT_i)
  for i in range(0, step - 1):
      k_1 = diffEq(RNAs[i], PROTs[i])
      # k_2 = diffEq(RNAs[i] + (deltaT/2), ys[i] + deltaT*(k_1/2))
      # k_3 = diffEq(ts[i] + (deltaT/2), ys[i] + deltaT*(k_2/2))
      # k_4 = diffEq(ts[i] + deltaT, ys[i] + deltaT*k_3)
      bRNAs.append(k_1[0]*deltaT + bRNAs[i])
      bPROTs.append(k_1[1]*deltaT + bPROTs[i])
  bRNAs = np.array(RNAs)
  bPROTs = np.array(PROTs)
  return [bRNAs, bPROTs, ts]

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
