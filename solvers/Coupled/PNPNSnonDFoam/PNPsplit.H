fvBlockMatrix<vector4> PNPsEqn(PNPs);

#include "psiEEqn.H"
#include "cEqn.H"
#include "couplingTermsPNPsplit.H"
//#include "calcF.H"

maxResidual = cmptMax(PNPsEqn.solve().initialResidual()); // residual = b - A*x_n

// Retrieve solution
PNPsEqn.retrieveSolution(0, cPlus.internalField());
PNPsEqn.retrieveSolution(1, cMinus.internalField());
PNPsEqn.retrieveSolution(2, psiI.internalField());
PNPsEqn.retrieveSolution(3, psiE.internalField());

cPlus.correctBoundaryConditions();
cMinus.correctBoundaryConditions();
psiI.correctBoundaryConditions();
psiE.correctBoundaryConditions();

#include "boundC.H"
#include "convergenceCheck.H"
