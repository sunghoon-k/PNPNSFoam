{
    Info<< "########################################### "  
        << "PNPEqn" 
        << " ###########################################" << nl << endl;
    
    fvBlockMatrix<vector3> PNPEqn(PNP);

    fvScalarMatrix psiEEqn
    (
        fvm::laplacian(psiE)
    );

    phiE = fvc::snGrad(psiE) * mesh.magSf();

    fvScalarMatrix cPlusEqn
    (
    - fvm::laplacian(DPlus_nonD, cPlus)
    - fvm::div(DPlus_nonD * phiE, cPlus)
    );

    fvScalarMatrix cMinusEqn
    (
    - fvm::laplacian(DMinus_nonD, cMinus)
    - fvm::div(DMinus_nonD * phiE, cMinus)
    );

    if(solveTransient)
    {
        cPlusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cPlus);
        cMinusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cMinus);
    }

    PNPEqn.insertEquation(0, psiEEqn);
    PNPEqn.insertEquation(1, cPlusEqn);
    PNPEqn.insertEquation(2, cMinusEqn);

    // Add off-diagonal coupling terms
    fvScalarMatrix Sp_fvScalarMatrix_cPlus(fvm::Sp(zPlus/sqr(DebL_nonD), cPlus));
    fvScalarMatrix Sp_fvScalarMatrix_cMinus(fvm::Sp(zMinus/sqr(DebL_nonD), cMinus));

    PNPEqn.insertEquationCoupling(0, 1, Sp_fvScalarMatrix_cPlus); // 4: psiE <- 5: cPlus
    PNPEqn.insertEquationCoupling(0, 2, Sp_fvScalarMatrix_cMinus); // 4: psiE <- 6: cMinus

    maxResidual_PNP = cmptMax(PNPEqn.solve().initialResidual()); // residual = b - A*x_n

    // Retrieve solution
    PNPEqn.retrieveSolution(0, psiE.internalField());
    PNPEqn.retrieveSolution(1, cPlus.internalField());
    PNPEqn.retrieveSolution(2, cMinus.internalField());

    psiE.correctBoundaryConditions();
    cPlus.correctBoundaryConditions();
    cMinus.correctBoundaryConditions();

    #include "boundC.H"
}
