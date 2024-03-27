fvBlockMatrix<vector3> PNPEqn(PNP);

#include "psiEEqn.H"
#include "cEqn.H"
#include "couplingTermsPNP.H"
//#include "calcF.H"

maxResidual = cmptMax(PNPEqn.solve().initialResidual()); // residual = b - A*x_n

// Retrieve solution
PNPEqn.retrieveSolution(0, cPlus.internalField());
PNPEqn.retrieveSolution(1, cMinus.internalField());
PNPEqn.retrieveSolution(2, psiE.internalField());
//PNPEqn.retrieveSolution(3, psiE.internalField());

cPlus.correctBoundaryConditions();
cMinus.correctBoundaryConditions();
psiE.correctBoundaryConditions();

#include "boundC.H"
#include "convergenceCheck.H"
