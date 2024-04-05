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

for Ratio in Ratios:
  for name in names:
    for Voltage in Voltages:

      for boundary in boundaries:
        addr_ = addr + f'/Transient/Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}'

        meshTitle = f"Diameter{Diameter}_Ratio{Ratio}_Dn{Dn-1}"
        os.system(f"xdotool key ctrl+shift+t")
        time.sleep(1)  
        os.system(f"xdotool type 'cd {addr_}/{meshTitle}_{boundary}_V{Voltage}'")
        time.sleep(0.5)        
        os.system(f"xdotool key Return")
        time.sleep(0.5)        
        os.system(f"xdotool type 'fe40'")
        time.sleep(0.5)        
        os.system(f"xdotool key Return")
        time.sleep(0.5)      
        os.system(f"xdotool type 'chmod 755 Allrun'")
        time.sleep(0.5)        
        os.system(f"xdotool key Return")
        time.sleep(0.5)          

        os.system(f"xdotool type './Allrun'") # './Allrun', 'PNPNSFoam_org_Ueof | tee log'
        time.sleep(0.5)        
        os.system(f"xdotool key Return")
        time.sleep(0.5)        

