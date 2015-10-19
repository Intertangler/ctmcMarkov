# Triangular Binding Simulation v8
from __future__ import division # this code enables the disply of 0.xxx numbers
# This is Kd1
import random

# Don't allow negative concentrations
def Add_Check(Amount, Change):
    Delta=Amount + Change
    if Delta >= 0.0:
        return Delta
    else:
        return 0.0
# Not needed here, but it's used to calculate the theoretical ka2
def Mono_To_Biv(ka1):
    Na=6.02*10**23
    Ab_Volume= 4/3*3.14*(14)**3 # nm3
    Ab_VolumeL= 0.5*Ab_Volume*10**-24 # 1 nm3 = 10^-24 L
    ka2=ka1*(1/(Ab_VolumeL*Na))
    return ka2
  
# A + B -> C, ka1 kd1 Monovalent, First Binding, Mono Site
# A + D -> E, ka1 kd1 Monovalent, First Binding, Biv Site
# E -> F, ka2 kd2 Bivalent, Second Binding, Biv Site

def Kd_Calculator_Tri(ka1,kd1,ka2, steps, ligands_per_sphere_1,ligands_per_sphere_2):
      Kd1=kd1/ka1
      ka2=ka2 # 1/Ms
      kd2=kd1 # 1/s
      Kd2=kd2/ka2
      N_A=1000*10**-9
      N_B=10*10**-9
      N_D=10*10**-9 #yo mamma!!!!!!!!
      N_Dall=N_D
      N_C=0.0
      N_E=0.0
      N_F=0.0
      l_s1=ligands_per_sphere_1
      l_s2=ligands_per_sphere_2
      for k in range(steps):
          N_Bdt = kd1*N_C - l_s1*ka1*N_A*N_B
          N_Cdt = -kd1*N_C + l_s1*ka1*N_A*N_B
          N_Ddt = kd1*N_E - l_s2*ka1*N_A*N_D
          N_Edt = -kd1*N_E + l_s2*ka1*N_A*N_D - ka2*(N_E/N_Dall)*(l_s2-1) + l_s2*kd2*N_F
          N_Fdt = ka2*(N_E/N_Dall)*(l_s2-1)-l_s2*kd2*N_F
          N_B = Add_Check(N_B, N_Bdt)
          N_C = Add_Check(N_C, N_Cdt)
          N_D = Add_Check(N_D, N_Ddt)
          N_E = Add_Check(N_E, N_Edt)
          N_F = Add_Check(N_F, N_Fdt)
      Kd_total = N_A*(l_s1*N_B+l_s2*N_D+(l_s2-1)*N_E)/(N_C+N_E+N_F)
      return [Kd1, Kd2, Kd_total]
      
def Kd_Calculator_Mono(ka1,kd1,steps):
      Kd1=kd1/ka1
      N_A=1000*10**-9
      N_B=10*10**-9
      N_C=0.0
      for k in range(steps):
          N_Bdt = kd1*N_C - ka1*N_A*N_B
          N_Cdt = -kd1*N_C + ka1*N_A*N_B
          N_B = Add_Check(N_B, N_Bdt)
          N_C = Add_Check(N_C, N_Cdt)
      Kd_total = N_A*(N_B)/(N_C)
      return [Kd1, Kd_total]

# Not needed here
def ka2_generator(initial_value,steps):
    ka2_list=[initial_value]
    for l in range(1,steps):
        ka2_temp=initial_value/(l*50)
        ka2_list.append(ka2_temp)
    return ka2_list



      
# Now let's try to simulate the real system
# A + B -> C, ka1 kd1 Monovalent, First Binding
# C + B -> E, ka2 kd2 Bivalent, Second Binding

#This is from the 1:1 model of a monovalent control
ka1=0.0001*1.74*10**7 # 1/Ms take steps of 0.1 milisecond
kd1=0.0001*5.08*10**-4 # 1/s take steps of 0.1 milisecond
ka2=4.87804878*10**-10 # 1/Ms take steps of 0.1 milisecond
results_list=[]

# This is the mono binding control
Kd_mono=Kd_Calculator_Mono(ka1,kd1,80000000)
temp=(Kd_mono,'Kd_Mono_Sim')
results_list.append(temp)

# this is the 2 monomeric antigen per structure control
results = Kd_Calculator_Tri(ka1,kd1,ka2,80000000,1,1)
temp=(results,'Kd_total_sim',ka2)
results_list.append(temp)

# this is the 1 bivalent antigen per structure control
results = Kd_Calculator_Tri(ka1,kd1,ka2,80000000,0,2)
temp=(results,'Kd_total_sim',ka2)
results_list.append(temp)

# this is the actual 1 monomeric 1 bivalent per structure simulation
results = Kd_Calculator_Tri(ka1,kd1,ka2,80000000,1,2)
temp=(results,'Kd_total_sim',ka2)
results_list.append(temp)
print results_list