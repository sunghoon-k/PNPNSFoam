import os
import subprocess
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

for Ratio in Ratios:
  for name in names:
    for Voltage in Voltages:

      str_head = """
shopt -s extglob
# RESULT=$(find ./ -name AR*)
# gmshToFoam $RESULT
rm -r !("0"|"constant"|"system"|"Allclean"|"Allrun")
# reset
PNPNSFoam_hoon | tee log.txt
"""
      for boundary in boundaries:
        meshTitle = f"Diameter{Diameter}_Ratio{Ratio}_Dn{Dn-1}"
        os.system(f"rm Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}/Allrun") # msh 파일 옮기기
        with open(f'Allrun', 'w') as f:
          f.write(str_head)
        os.system(f"mv Allrun Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}/") # msh 파일 옮기기
