import os
import subprocess
import numpy as np 

Voltages = [Volt for Volt in range(20, 111, 30)] # [Volt for Volt in range(21, 41)]
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
    for Voltage in Voltages:

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
    location    "constant";
    object      transportProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
transportModel Newtonian;
nu              nu [0 2 -1 0 0 0 0] 1;

phiDyMBoundary
{
	name	top;
"""
      str_head = str_head + f"""
	phiStart	{Voltage};
	phiEnd		{Voltage};
"""
      str_head = str_head + """
	phiInterval	1;
}

Ions
{
    ionsList(cPlus,cMinus); 

	cPlus
	{
		variableName cPlus; 		
		z z [0 0 0 0 0 0 0] 1; // valence
		D D [0 -3 -1 0 0 0 0] 1.33e-9; // diffusivity
	}

	cMinus
	{
		variableName cMinus; 		
		z z [0 0 0 0 0 0 0] -1; // valence
		D D [0 -3 -1 0 0 0 0] 2.03e-9; // diffusivity
 	}
}

Electrolyte
{
  	Cation C1;
	Anion C2; 

	Potential Phi;
	Velocity U;
  	Pressure p;

"""
      str_body = str_head + f"""    l0              l0 [0 1 0 0 0 0 0] {Height}e-6; // length-scale 

"""
      str_tail = str_body + """
    c0              c0 [0 -3 0 0 1 0 0] 10; // Bulk Concentration   
    T               T  [0 0 0 1 0 0 0] 300; // Absolute temperature    
    e               e  [0 0 1 0 0 1 0] 1.60217646e-19; // elementary charge
    kB              kB [1 2 -2 -1 0 0 0] 1.38065e-23; // Bolt-Zmann constant
    F               F  [0 0 1 0 -1 1 0] 96485.3415; // Faraday constant
    NA              NA [0 0 0 0 -1 0 0] 6.0221415e23; // Avogadro's number 
    eps0          eps0 [0 -1 1 0 -1 1 0] 8.8541878e-12; // vacumm permittivity
    epsr          epsr [0 0 0 0 0 0 0] 80; // dielectric constant 
    mu              mu [1 -1 -1 0 0 0 0] 1e-3; // dynamic viscosity
    rho            rho [1 -3 0 0 0 0 0] 1000; // density
    nu              nu [0 2 -1 0 0 0 0] 1; // kinematic viscosity 
    wsc            wsc [0 -2 1 0 0 1 0] -2e-3 ; // Wall's Surface Charge   
}

// ************************************************************************* //    
"""
      meshTitle = f"Diameter{Diameter}_Ratio{Ratio}_Dn{Dn-1}"

      for boundary in boundaries:
        with open(f'transportProperties', 'w') as f:
          f.write(str_tail)
        os.system(f"mv transportProperties Diameter{Diameter}/Ratio{Ratio}_Dn{Dn-1}/V{Voltage}/{meshTitle}_{boundary}_V{Voltage}/constant/")
