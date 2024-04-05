import os
import subprocess
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

    str_head = """
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  5                                     |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
application     N-S-PNPDImlessFoam;
startFrom       latestTime; // startTime;
startTime       0;
stopAt          endTime;
endTime         100; // 1;
deltaT          1;
writeControl    timeStep;
writeInterval   1;
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
libs ("libfixedIonicFlux.so" "libmodifiedFixedValue.so" "libslipSmoluchowski.so" "libzeroIonicFlux_nonD.so");
nNewtonIter	50; //15;
nNewtonIterPNP	150;
nPNPNSIter	500;
C1C2convergence	1e-3; //1e-3;
PNPNStransient no;
solveNS no; // yes;
pseudoTransient no;
splitpsi no;
saveVtkFiles no;
saveDesaltPhysics no;
debugMode no;
// ************************************************************************* //
"""
    meshTitle = f"Diameter{Diameter}_Ratio{Ratio}_Dn{Dn-1}"

    for boundary in boundaries:
      with open(f'controlDict', 'w') as f:
        f.write(str_head)
      os.system(f"mv controlDict Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/{meshTitle}_{boundary}/system/")
