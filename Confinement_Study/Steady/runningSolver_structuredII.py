import os
import subprocess
import time
import numpy as np 

boundaries = ['slip', 'noSlip'] # ['slip', 'Smoluchowski']
Ratios  = [2, 4, 8, 16] # [1, 2, 4, 8, 16] # 종횡비
Diameter = 20
Voltage = 150
names = ['Diameter'] # ['Diameter', 'Square', 'Square_recombine']

dl = 1/3 # 1/3 um
dr = np.sqrt(4/np.sqrt(3)) * dl
Dn = int(Diameter*np.pi/dr/4)
Dn = int(Dn/2) + 1 # coarse: 기존의 1/2 만큼

for Ratio in Ratios:
  for name in names:
    Height = Diameter / Ratio
    
    for boundary in boundaries:

      addr = f'/media/sh/drive3/24_04/Steady/Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}'
  
      meshTitle = f"Diameter{Diameter}_Ratio{Ratio}_Dn{Dn-1}"
      os.system(f"xdotool key ctrl+shift+t")
      time.sleep(1)  
      os.system(f"xdotool type 'cd {addr}/{meshTitle}_{boundary}'")
      time.sleep(0.5)        
      os.system(f"xdotool key Return")
      time.sleep(0.5)        
      os.system(f"xdotool type 'fe40'")
      time.sleep(0.5)        
      os.system(f"xdotool key Return")
      time.sleep(0.5)        
      os.system(f"xdotool type './Allrun'") # './Allrun', 'PNPNSFoam_org_Ueof > log &'
      time.sleep(0.5)        
      os.system(f"xdotool key Return")
      time.sleep(0.5)        

