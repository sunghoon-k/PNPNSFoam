fvBlockMatrix<vector8> PNPNSsEqn(PNPNSs);

p.storePrevIter();

#include "cEqn.H"
#include "psiEEqn.H"
#include "UEqn.H"
#include "pEqn.H"
#include "couplingTermsPNPNSsplit.H"

Info << "before solve" << endl;

maxResidual = cmptMax(PNPNSsEqn.solve().initialResidual()); // residual = b - A*x_n

// Retrieve solution
PNPNSsEqn.retrieveSolution(0, cPlus.internalField());
PNPNSsEqn.retrieveSolution(1, cMinus.internalField());
PNPNSsEqn.retrieveSolution(2, psiI.internalField());
PNPNSsEqn.retrieveSolution(3, psiE.internalField());
PNPNSsEqn.retrieveSolution(4, U.internalField());
PNPNSsEqn.retrieveSolution(7, p.internalField());

cPlus.correctBoundaryConditions();
cMinus.correctBoundaryConditions();
psiI.correctBoundaryConditions();
psiE.correctBoundaryConditions();
U.correctBoundaryConditions();
p.correctBoundaryConditions();

phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;
#include "continuityErrs.H"
p.relax();

#include "boundC.H"
#include "convergenceCheck.H"
