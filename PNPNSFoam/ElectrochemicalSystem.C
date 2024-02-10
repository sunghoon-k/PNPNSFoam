#include "ElectrochemicalSystem.H"
#include "wordList.H"
#include "IFstream.H"
#include <iostream>
#include <iomanip>
namespace Foam
{
	defineTypeNameAndDebug(ElectrochemicalSystem, 0);
}

// independent constant

/*
e(dimensionedScalar("e",		dimensionSet( 0, 0, 1, 0, 0, 1, 0 ),	1.6021766208e-19)),
kB(dimensionedScalar("kB",		dimensionSet( 1, 2, -2, -1, 0, 0, 0 ),	1.38064852e-23)),
eps0(dimensionedScalar("eps0",	dimensionSet(-1, -3, 4, 0, 0, 2, 0),	8.8541878176e-12)),
NA(dimensionedScalar("NA",		dimensionSet(0, 0, 0, 0, -1, 0, 0),		6.022140857e+23)),
F(dimensionedScalar("F",		(e * NA).dimensions(),					(e * NA).value())),
*/


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
objects(mesh_, mesh_.time().timeName()),

// independent constant
e("e", dimensionSet( 0, 0, 1, 0, 0, 1, 0 ),	1.6021766208e-19),
kB("kB", dimensionSet( 1, 2, -2, -1, 0, 0, 0 ),	1.38064852e-23),
eps0("eps0", dimensionSet(-1, -3, 4, 0, 0, 2, 0),	8.8541878176e-12),
NA("NA", dimensionSet(0, 0, 0, 0, -1, 0, 0),		6.022140857e+23),
F("F", (e * NA).dimensions(),	(e * NA).value()),

// dependent constant
T(dict.subDict("Electrolyte").lookup("T")),
l0(dict.subDict("Electrolyte").lookup("l0")),
c0(dict.subDict("Electrolyte").lookup("C0")),
epsr(dict.subDict("Electrolyte").lookup("epsr")),
// wsc(dict.subDict("Electrolyte").lookup("wsc")),
rho(dict.subDict("Electrolyte").lookup("rho")),
// mu(dict.subDict("Electrolyte").lookup("mu")),
nu(dict.subDict("Electrolyte").lookup("nu")),
zPlus(dict.subDict("Ions").subDict("cPlus").lookup("zPlus")),
DPlus(dict.subDict("Ions").subDict("cPlus").lookup("DPlus")),
zMinus(dict.subDict("Ions").subDict("cMinus").lookup("zMinus")),
DMinus(dict.subDict("Ions").subDict("cMinus").lookup("DMinus"))

{
  /*
  // dependent constant
  dimensionedScalar Foam::ElectrochemicalSystem::T(dict.subDict("Electrolyte").lookup("T"));
  dimensionedScalar Foam::ElectrochemicalSystem::l0(dict.subDict("Electrolyte").lookup("l0"));
  dimensionedScalar Foam::ElectrochemicalSystem::c0(dict.subDict("Electrolyte").lookup("C0"));
  dimensionedScalar Foam::ElectrochemicalSystem::epsr(dict.subDict("Electrolyte").lookup("epsr"));
  // wsc(dict.subDict("Electrolyte").lookup("wsc")),
  dimensionedScalar Foam::ElectrochemicalSystem::rho(dict.subDict("Electrolyte").lookup("rho"));
  // mu(dict.subDict("Electrolyte").lookup("mu")),
  dimensionedScalar Foam::ElectrochemicalSystem::nu(dict.subDict("Electrolyte").lookup("nu"));
  dimensionedScalar Foam::ElectrochemicalSystem::zPlus(dict.subDict("Ions").subDict("cPlus").lookup("zPlus"));
  dimensionedScalar Foam::ElectrochemicalSystem::DPlus(dict.subDict("Ions").subDict("cPlus").lookup("DPlus"));
  dimensionedScalar Foam::ElectrochemicalSystem::zMinus(dict.subDict("Ions").subDict("cMinus").lookup("zMinus"));
  dimensionedScalar Foam::ElectrochemicalSystem::DMinus(dict.subDict("Ions").subDict("cMinus").lookup("DMinus"));

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


