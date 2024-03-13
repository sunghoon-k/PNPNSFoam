#ifndef ElectrochemicalSystem_H
#define ElectrochemicalSystem_H

#include "fvCFD.H"
#include "vtkMesh.H"
#include "cmtVTU.H"
#include "internalWriter.H"
#include "wordReList.H"
#include "IOobjectList.H"
#include "readFields.H"
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
	word Mem1Name,Mem2Name;///* Sang add
	word U_,C1_,C2_,Phi_,P_;///* Sang add
	word desaltType;
	wordList desalinationInletsList;
	wordList desalinationOutletsList;
  bool writeHeader_;
	bool solveNS_;
	IOobjectList objects;
public:
  dimensionedScalar D1, D2;
  dimensionedScalar z1, z2;
  dimensionedScalar l0, C0, T, e, kB, F, NA, eps0, epsr, wsc, rho, mu;
  dimensionedScalar nu;
  scalar DimlD1, DimlD2, Pe;

	//- Private member functions
    template<class GeoField> void print
    (
    	const char* msg,
    	Ostream& os,
    	const PtrList<GeoField>& flds
    );
	void print
	(
		Ostream& os,
		const wordList& flds
	);

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
  void changePhi();
  void saveDataVTK(const fvMesh& mesh,int timestep=-1);
  void saveDataVTK(string filnename, PtrList<volScalarField> vsf);
  void saveDesaltPhysics();
  void saveIV();
  void CalPowerDissipation(string filenamehead);
  scalar calcIonicFluxJ(word patchName, word fieldName, scalar DimlD, dimensionedScalar z);
  void calcC(word fieldName, scalar& qc1, scalar& qc2, wordList physicalpatchName);
  scalar I1(word patchName);
  scalar I2(word patchName);
  scalar I(wordList physicalpatchName,word C12);
  scalar CurrentSum();
  scalar massflux(label patchID,label faceID);
  scalar CurrentPerIteration(scalar i,bool writeHeader);
  void CurrentPerTimestep(scalar timestep);
  void writeIonicFluxBoundary(word fieldName,scalar D,scalar z,scalar phiInstant,word patchName);

}; // End class

} // End namespace Foam

#endif
