#include "ElectrochemicalSystem.H"
#include "wordList.H"
#include "IFstream.H"
#include <iostream>
#include <iomanip>
namespace Foam
{
	defineTypeNameAndDebug(ElectrochemicalSystem, 0);
}

// constructor
Foam::ElectrochemicalSystem::ElectrochemicalSystem
(
	const fvMesh& mesh,
	const dictionary& dict
):
mesh_(mesh),
electricdict_(dict),
phiStart_(readScalar(dict.subDict("phiDyMBoundary").lookup("phiStart"))),
phiEnd_(readScalar(dict.subDict("phiDyMBoundary").lookup("phiEnd"))),
phiInterval_(readScalar(dict.subDict("phiDyMBoundary").lookup("phiInterval"))),
phiInstant_(phiStart_),
objects(mesh_, mesh_.time().timeName())
/*
l0(dict.subDict("Electrolyte").lookup("l0")),
c0(dict.subDict("Electrolyte").lookup("C0")),
T_const(dict.subDict("Electrolyte").lookup("T_const")),
e_const(dict.subDict("Electrolyte").lookup("e_const")),
kB_const(dict.subDict("Electrolyte").lookup("kB_const")),
F_const(dict.subDict("Electrolyte").lookup("F_const")),
NA_const(dict.subDict("Electrolyte").lookup("NA_const")), // Avogadro
eps0(dict.subDict("Electrolyte").lookup("eps0")),
epsr(dict.subDict("Electrolyte").lookup("epsr")),
// wsc(dict.subDict("Electrolyte").lookup("wsc")),
rho(dict.subDict("Electrolyte").lookup("rho")),
// mu(dict.subDict("Electrolyte").lookup("mu")),
nu(dict.subDict("Electrolyte").lookup("nu")),
zPlus(dict.subDict("Ions").subDict("cPlus").lookup("zPlus")),
DPlus(dict.subDict("Ions").subDict("cPlus").lookup("DPlus")),
zMinus(dict.subDict("Ions").subDict("cMinus").lookup("zMinus")),
DMinus(dict.subDict("Ions").subDict("cMinus").lookup("DMinus"))
*/
{
  /*
  scalar D0 = 0.5*(DPlus.value() + DMinus.value());
  DimlDPlus = DPlus.value()/D0;
  DimlDMinus = DMinus.value()/D0;

  scalar psiE0 = (kB_const.value()*T_const.value())/(zPlus.value()*e_const.value()); // Thermal Voltage
  scalar U0 = epsr.value()*eps0.value()*sqr(Phi0)/(mu.value()*l0.value());
  scalar Pe = U0*l0.value()/D0; // Peclet #
  scalar P0 = nu.value();
  scalar t0 = sqrt(l0)/D0; 

  mu = nu/rho; 
  Info<< "ElectricProperties" << endl;
  Info << "****************************************************************" << endl;
  Info << "length-scale: " << l0.value() << endl;
  Info << "velocity-scale: " << U0.value() << endl;
  Info << "time-scale: " << t0.value() << endl;
  Info << "Thermal voltage: " << psiE0.value() << endl;
  Info << "Bulk Concentration :" << c0.value() << endl;
  Info << "Absolute temperature: " << T_const.value() << endl;
  Info << "elementary charge: " << e_const.value() << endl;
  Info << "Bolt-Zmann constant: " << kB_const.value() << endl;
  Info << "Faraday constant: " << F_const.value() << endl;
  Info << "Avogadro's number: " << NA_const.value() << endl;
  Info << "Vacumm permittivity: " << eps0.value() << endl;
  Info << "Dielectric constant: " << epsr.value() << endl;
//  Info << "Wall Surface Charge: " << wsc.value() << endl;
  Info << "Density: " << rho.value() << endl;
  Info << "Dynamics viscosity: " << mu.value() << endl; // mu = nu/rho [m^2/s]
  Info << "Kinematic viscosity: " << nu.value() << endl; // nu [N s/m^2]
  Info << "Diffusivity DPlus: " << DPlus.value() << endl;
  Info << "Diffusivity DMinus: " << DMinus.value() << endl;
  */
}


