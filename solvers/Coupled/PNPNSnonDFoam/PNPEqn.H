{
    fvBlockMatrix<vector3> PNPEqn(PNP);

    fvScalarMatrix psiEEqn
    (
        fvm::laplacian(psiE)
    );

    phiE = fvc::snGrad(psiE) * mesh.magSf();

    fvScalarMatrix cPlusEqn
    (
    - fvm::laplacian(DPlus, cPlus)
    - fvm::div(mobilityPlus * phiE, cPlus)
    );

    fvScalarMatrix cMinusEqn
    (
    - fvm::laplacian(DMinus, cMinus)
    - fvm::div(mobilityMinus * phiE, cMinus)
    );

    if(solveTransient)
    {
        cPlusEqn += fvm::ddt(cPlus);
        cMinusEqn += fvm::ddt(cMinus);
    }

    PNPEqn.insertEquation(0, psiEEqn);
    PNPEqn.insertEquation(1, cPlusEqn);
    PNPEqn.insertEquation(2, cMinusEqn);

    // Add off-diagonal coupling terms
    fvScalarMatrix Sp_fvScalarMatrix_cPlus(fvm::Sp((e*zPlus)/(eps0*epsr), cPlus));
    fvScalarMatrix Sp_fvScalarMatrix_cMinus(fvm::Sp((e*zMinus)/(eps0*epsr), cMinus));

    PNPEqn.insertEquationCoupling(0, 1, Sp_fvScalarMatrix_cPlus); // 4: psiE <- 5: cPlus
    PNPEqn.insertEquationCoupling(0, 2, Sp_fvScalarMatrix_cMinus); // 4: psiE <- 6: cMinus

    PNPEqn.solve(); // residual = b - A*x_n

    // Retrieve solution
    PNPEqn.retrieveSolution(0, psiE.internalField());
    PNPEqn.retrieveSolution(1, cPlus.internalField());
    PNPEqn.retrieveSolution(2, cMinus.internalField());

    psiE.correctBoundaryConditions();
    cPlus.correctBoundaryConditions();
    cMinus.correctBoundaryConditions();
}
