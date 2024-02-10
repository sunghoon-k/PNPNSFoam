#ifndef ElectrochemicalSystem_H
#define ElectrochemicalSystem_H

#include "fvCFD.H"
//#include "vtkMesh.H"
//#include "cmtVTU.H"
//#include "internalWriter.H"
#include "wordReList.H"
#include "IOobjectList.H"
// #include "readFields.H"
#include "pointMesh.H"
#include <stdlib.h>
#include <stdio.h>
namespace Foam
{

class ElectrochemicalSystem
{
	const fvMesh& mesh_;
	const dictionary& electricdict_;
	scalar phiStart_;
	scalar phiEnd_;
	scalar phiInterval_;
	scalar phiInstant_;
	IOobjectList objects;
/*
public:
	dimensionedScalar DPlus, DMinus;
	dimensionedScalar zPlus, zMinus;
	dimensionedScalar l0, c0, T_const, e_const, kB_const, F_const, NA_const, eps0, epsr, rho, mu, nu; // wsc
	scalar DimlD1, DimlD2, Pe; // dimensionless constant
*/

public:
	TypeName("ElectrochemicalSystem");

	ElectrochemicalSystem
	(
		const fvMesh& mesh,
		const dictionary& dict
	);

	virtual ~ElectrochemicalSystem()
	{}

	inline bool phiRun()
	{
		phiStart_ = readScalar(electricdict_.subDict("phiDyMBoundary").lookup("phiStart"));
		phiEnd_ = readScalar(electricdict_.subDict("phiDyMBoundary").lookup("phiEnd"));
		phiInterval_ = readScalar(electricdict_.subDict("phiDyMBoundary").lookup("phiInterval"));
		if(mag(phiInstant_) > mag(phiEnd_))
			return false;
		else
			return true;
	};

	inline scalar phiInstant()
	{
		return phiInstant_;
	};
	
	void set_phiInstant(scalar phi)
	{
		phiInstant_ = phi;
	};

	void changePhi()
	{
		phiInstant_ += phiInterval_;
	}


}; // End class

} // End namespace Foam

#endif
