p.storePrevIter();

fvBlockMatrix<vector7> PNPNSEqn(PNPNS);

netCharge = (zPlus*cPlus + zMinus*cMinus);

#include "UEqn.H"
#include "pEqn.H"
#include "psiEEqn.H"
#include "cEqn.H"
#include "couplingTerms.H"

maxResidual = cmptMax(PNPNSEqn.solve().initialResidual()); // residual = b - A*x_n

// Retrieve solution
PNPNSEqn.retrieveSolution(0, U.internalField());
PNPNSEqn.retrieveSolution(3, p.internalField());
PNPNSEqn.retrieveSolution(4, psiE.internalField());
PNPNSEqn.retrieveSolution(5, cPlus.internalField());
PNPNSEqn.retrieveSolution(6, cMinus.internalField());

U.correctBoundaryConditions();
p.correctBoundaryConditions();
psiE.correctBoundaryConditions();
cPlus.correctBoundaryConditions();
cMinus.correctBoundaryConditions();

phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;

#include "continuityErrs.H"
#include "boundC.H"

p.relax();

#include "convergenceCheck.H"