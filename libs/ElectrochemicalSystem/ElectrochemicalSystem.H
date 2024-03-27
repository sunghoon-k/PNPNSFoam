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

public:
	// independent constant
	static dimensionedScalar e, kB, eps0, NA, F;

	// dependent constant
	dimensionedScalar T, l0, c0, epsr, rho, mu, zPlus, DPlus, zMinus, DMinus, psiE0;
	// scalar DimlD1, DimlD2, Pe; // dimensionless constant

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

	inline scalar phiStart()
	{
		return phiStart_;
	};

	inline scalar phiInterval()
	{
		return phiInterval_;
	};

	inline scalar phiEnd()
	{
		return phiEnd_;
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
