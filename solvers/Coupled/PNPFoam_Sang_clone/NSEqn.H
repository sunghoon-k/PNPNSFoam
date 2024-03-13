#include "readBlockSolverControls.H"

p.storePrevIter();
fvBlockMatrix<vector4> UpEqn(Up);

#include "UEqn.H"
#include "pEqn.H"
#include "couplingTerms.H"

// Solve the block matrix
UpEqn.solve();

// Retrieve solution
UpEqn.retrieveSolution(0, U.internalField());
UpEqn.retrieveSolution(3, p.internalField());

U.correctBoundaryConditions();
p.correctBoundaryConditions();

phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;
#include "continuityErrs.H"
p.relax();
turbulence->correct();
