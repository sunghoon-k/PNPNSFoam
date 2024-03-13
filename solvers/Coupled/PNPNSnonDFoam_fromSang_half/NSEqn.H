bodyForce = - netCharge * fvc::grad(psiE) * sqrDebL_nonD/2;
nNewtonIter = 1;
for(int NewtonIter = 0; NewtonIter < nNewtonIter; NewtonIter++)
{
    p.storePrevIter();

    fvBlockMatrix<vector4> NSEqn(NS);

    #include "UEqn.H"
    #include "pEqn.H"

    BlockLduSystem<vector, vector> pInU(fvm::grad(p));
    NSEqn.insertBlockCoupling(0, 3, pInU, true); // 0: U <- 3: p

    BlockLduSystem<vector, scalar> UInp(fvm::UDiv(U));
    NSEqn.insertBlockCoupling(3, 0, UInp, false); // 3: p <- 0: U

    maxResidual = cmptMax(NSEqn.solve().initialResidual()); // residual = b - A*x_n

    // Retrieve solution
    NSEqn.retrieveSolution(0, U.internalField());
    NSEqn.retrieveSolution(3, p.internalField());

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;

    #include "continuityErrs.H"
    // #include "boundPU.H"

    p.relax();

    #include "convergenceCheck.H"
}
