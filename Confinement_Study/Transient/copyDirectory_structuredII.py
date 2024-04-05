import os
import subprocess
import time
import numpy as np

Voltages = [Volt for Volt in range(20, 111, 30)] # [Volt for Volt in range(21, 41)]
boundaries = ['slip', 'noSlip'] # ['slip', 'Smoluchowski']
Ratios = [2, 4, 8, 16] # [1, 2, 8, 16]
Diameter = 20
addr = '/media/sh/drive3/24_04'
names = ['Diameter'] # ['Diameter', 'Square', 'Square_recombine']

dl = 1/3 # 1/3 um
dr = np.sqrt(4/np.sqrt(3)) * dl
Dn = int(Diameter*np.pi/dr/4)
Dn = int(Dn/2) + 1 # coarse: 기존의 1/2 만큼

os.system(f"mkdir Diameter{Diameter}")
for Ratio in Ratios:
  Height = Diameter / Ratio
  for name in names:
    os.system(f"mkdir Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}")

    for Voltage in Voltages:
      os.system(f"mkdir Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}")
  
      for boundary in boundaries:
        meshTitle = f"Diameter{Diameter}_Ratio{Ratio}_Dn{Dn-1}"
        print(f"{addr}/Transient/Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}")
        os.system(f"cp -r {addr}/Steady/Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/{meshTitle}_{boundary} {addr}/Transient/Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}")
        os.system(f"rm -r {addr}/Transient/Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}/0")
        os.system(f"mv {addr}/Transient/Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}/{Voltage+1} {addr}/Transient/Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}/0") # Voltage+1 or 1

