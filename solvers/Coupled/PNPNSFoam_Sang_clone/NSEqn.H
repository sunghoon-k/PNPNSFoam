volScalarField netCharge = C1*z1.value() + C2*z2.value();
volVectorField bodyForce = - netCharge * fvc::grad(Phi) / (Epsilon2*2);
maxResidual = 10;

for(int NewtonIter = 0; NewtonIter < nNewtonIteration; NewtonIter++)
{
    #include "readBlockSolverControls.H"

    p.storePrevIter();
    fvBlockMatrix<vector4> UpEqn(Up);

    #include "UEqn.H"
    #include "pEqn.H"
    #include "couplingTerms.H"

    // Solve the block matrix
    //UpEqn.solve();
    maxResidual = cmptMax(UpEqn.solve().initialResidual()); // residual = b - A*x_n

    // Retrieve solution
    UpEqn.retrieveSolution(0, U.internalField());
    UpEqn.retrieveSolution(3, p.internalField());

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;
    #include "continuityErrs.H"
    p.relax();
    turbulence->correct();

    #include "convergenceCheck.H"
}