import os
import subprocess
import time

import numpy as np
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

def perturb(value):
    perturbation = np.random.uniform(-0.01, 0.01)  # generates a random number between -0.01 and 0.01
    return value * (1 + perturbation)  # applies the perturbation to the given value

def perturb_field(field):
    return [perturb(value) for value in field]

ions = ['C1', 'C2']

Voltages = [Volt for Volt in range(20, 151, 10)] # [Volt for Volt in range(21, 41)]
boundaries = ['Smoluchowski'] # ['slip', 'Smoluchowski']
Ratios = [4] # [1, 2, 8, 16]
Height = 5
addr = '/media/sh/drive3/October2ndWeek'
names = ['Diameter'] # ['Diameter', 'Square', 'Square_recombine']

for Ratio in Ratios:
  dl = 1/3 # 1/3 um
  dr = np.sqrt(4/np.sqrt(3)) * dl
  Diameter = Height * Ratio
  Dn = int(Diameter*np.pi/dr/4)
  Dn = int(Dn/2) + 1 # coarse: 기존의 1/2 만큼

  for name in names:
    for Voltage in Voltages:
      for boundary in boundaries:
        meshTitle = f"Height{Height}_Ratio{Ratio}_Dn{Dn-1}"
        print(f"{addr}/Transient_c01/Height{Height}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}")

        for ion in ions:
          filename = f"{addr}/Transient_c01/Height{Height}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}/0/{ion}"
          foam_file = ParsedParameterFile(filename)
            
          # Let's say you want to perturb the 'value' field in the foam_file
          original_field = foam_file['internalField']  # replace 'value' with the name of the field you want to perturb
          perturbed_field = perturb_field(original_field)
    
          # Update the value in the file
          foam_file['internalField'] = perturbed_field  # replace 'value' with the name of the field you want to perturb
      
          # Write the file back in the original format
          foam_file.writeFile()

          # Read the file contents
          with open(filename, 'r') as file:
            lines = file.readlines()

          # Add a comment to the 15th line
          lines[15-1] = 'internalField   nonuniform List<scalar>\n'  # Assuming line numbering starts from 1

          # Write the file contents back to the file
          with open(filename, 'w') as file:
            file.writelines(lines)
  
