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
    os.system(f"mkdir {addr}/Result")
    os.system(f"mkdir {addr}/Result/Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}")

    for Voltage in Voltages:
      os.system(f"mkdir {addr}/Result/Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}")

      for boundary in boundaries:
      # os.system(f"mkdir {addr}/Result/Diameter{Diameter}/Ratio{Ratio}/boundary{boundary}")          

        meshTitle = f"Diameter{Diameter}_Ratio{Ratio}_Dn{Dn-1}"
        os.system(f"touch Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}/{meshTitle}_{boundary}_V{Voltage}.foam")
        time.sleep(1)


'''
for Height in Heights:
  os.system(f"mkdir {addr}/result/Height{Height}")
  for Ratio in Ratios:
    os.system(f"mkdir {addr}/result/Height{Height}/Ratio{Ratio}")
    for Ratio1 in Ratio1s:
      os.system(f"mkdir {addr}/result/Height{Height}/Ratio{Ratio}/mesh{Ratio1}")
      for Voltage in Voltages:
        for boundary in boundaries:
          os.system(f"mkdir {addr}/result/Height{Height}/Ratio{Ratio}/mesh{Ratio1}/boundary{boundary}")          
          D = Height * Ratio
          meshTitle = f"D{D}H{Height}_Ratio{str(Ratio1).replace('.', '')}"

          os.system(f"touch Height{Height}/Ratio{Ratio}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}/{meshTitle}_{boundary}_V{Voltage}.foam")
          time.sleep(1)

          os.system(f"xdotool key ctrl+shift+t")
          time.sleep(1)
          os.system(f"xdotool type 'cd {addr}/Height{Height}/Ratio{Ratio}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}'")
          os.system(f"xdotool key Return")
          os.system(f"xdotool type 'touch {meshTitle}_{boundary}_V{Voltage}.foam'")
          os.system(f"xdotool key Return")
'''
