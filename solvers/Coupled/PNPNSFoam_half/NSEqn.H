Info<< "########################################### "
    << "NSEqn"
    << " ###########################################" << nl << endl;            

p.storePrevIter();

fvBlockMatrix<vector4> NSEqn(NS);

netCharge = e*(zPlus*cPlus + zMinus*cMinus);

#include "UEqn.H"
#include "pEqn.H"
#include "couplingTerms.H"

maxResidual_NS = cmptMax(NSEqn.solve().initialResidual()); // residual = b - A*x_n

// Retrieve solution
NSEqn.retrieveSolution(0, U.internalField());
NSEqn.retrieveSolution(3, p.internalField());

U.correctBoundaryConditions();
p.correctBoundaryConditions();

phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;

#include "continuityErrs.H"
#include "boundPU.H"

p.relax();

//#include "convergenceCheck.H"