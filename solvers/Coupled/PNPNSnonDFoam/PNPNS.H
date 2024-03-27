fvBlockMatrix<vector7> PNPNSEqn(PNPNS);

p.storePrevIter();

#include "cEqn.H"
#include "psiEEqn.H"
#include "UEqn.H"
#include "pEqn.H"
#include "couplingTermsPNPNS.H"

maxResidual = cmptMax(PNPNSEqn.solve().initialResidual()); // residual = b - A*x_n

// Retrieve solution
PNPNSEqn.retrieveSolution(0, cPlus.internalField());
PNPNSEqn.retrieveSolution(1, cMinus.internalField());
PNPNSEqn.retrieveSolution(2, psiE.internalField());
PNPNSEqn.retrieveSolution(3, U.internalField());
PNPNSEqn.retrieveSolution(6, p.internalField());

cPlus.correctBoundaryConditions();
cMinus.correctBoundaryConditions();
psiE.correctBoundaryConditions();
U.correctBoundaryConditions();
p.correctBoundaryConditions();

phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;
#include "continuityErrs.H"
p.relax();

#include "boundC.H"
#include "convergenceCheck.H"
