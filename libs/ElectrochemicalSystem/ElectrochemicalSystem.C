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

dimensionedScalar ElectrochemicalSystem::e(dimensionedScalar("e",		dimensionSet( 0, 0, 1, 0, 0, 1, 0 ),	1.6021766208e-19));
dimensionedScalar ElectrochemicalSystem::kB(dimensionedScalar("kB",		dimensionSet( 1, 2, -2, -1, 0, 0, 0 ),	1.38064852e-23));
dimensionedScalar ElectrochemicalSystem::eps0(dimensionedScalar("eps0",	dimensionSet(-1, -3, 4, 0, 0, 2, 0),	8.8541878176e-12));
dimensionedScalar ElectrochemicalSystem::NA(dimensionedScalar("NA",		dimensionSet(0, 0, 0, 0, -1, 0, 0),		6.022140857e+23));
dimensionedScalar ElectrochemicalSystem::F(dimensionedScalar("F",		(e * NA).dimensions(),					(e * NA).value()));

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

T(dict.subDict("Electrolyte").lookup("T")),
l0(dict.subDict("Electrolyte").lookup("l0")),
c0(dict.subDict("Electrolyte").lookup("c0")),
epsr(dict.subDict("Electrolyte").lookup("epsr")),
// wsc(dict.subDict("Electrolyte").lookup("wsc")),
rho(dict.subDict("Electrolyte").lookup("rho")),
// mu(dict.subDict("Electrolyte").lookup("mu")),
mu(dict.subDict("Electrolyte").lookup("mu")),
zPlus(dict.subDict("Ions").subDict("cPlus").lookup("z")),
DPlus(dict.subDict("Ions").subDict("cPlus").lookup("D")),
zMinus(dict.subDict("Ions").subDict("cMinus").lookup("z")),
DMinus(dict.subDict("Ions").subDict("cMinus").lookup("D")),
psiE0(kB*T/e)
{
}


