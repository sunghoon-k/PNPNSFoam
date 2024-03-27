#include "readBlockSolverControls.H"
fvBlockMatrix<vector4> UpEqn(Up);

for(int NewtonIter = 0; NewtonIter < nNewtonIterNS; NewtonIter++)
{
    p.storePrevIter();

    #include "UEqn.H"
    #include "pEqn.H"
    #include "couplingTermsNS.H"

    // Solve the block matrix
    maxResidual_NS = cmptMax(UpEqn.solve().initialResidual());

    // Retrieve solution
    UpEqn.retrieveSolution(0, U.internalField());
    UpEqn.retrieveSolution(3, p.internalField());

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;
    #include "continuityErrs.H"
    Info << "maxResidual in NS: " << maxResidual_NS << endl;
    if (maxResidual_NS < convergenceCriterionNS)
    {
        Info << "reached criterion!" << endl; 
        break; 
    }

    p.relax();
}
