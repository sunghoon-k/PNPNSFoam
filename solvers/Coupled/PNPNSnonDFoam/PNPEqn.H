{
    bool bounded = false;

    fvBlockMatrix<vector3> PNPEqn(PNP);

    fvScalarMatrix psiEEqn
    (
        fvm::laplacian(psiE)
        //==
        //- netCharge/(eps0*epsr)
    );
    PNPEqn.insertEquation(0, psiEEqn);

    phiE = fvc::snGrad(psiE) * mesh.magSf();

    DPlusphiE = DPlus_nonD*zPlus*phiE;
    DMinusphiE = DMinus_nonD*zMinus*phiE;

    fvScalarMatrix cPlusEqn
    (
    - fvm::laplacian(DPlus_nonD, cPlus)
    - fvm::div(DPlusphiE, cPlus)
    //- fvm::div(DPlus_nonD*zPlus*phiE, cPlus) // 
    //- fvm::div(mobilityPlus * phiE, cPlus)
    //- fvc::laplacian(mobilityPlus * cPlus, psiE)
    );

    fvScalarMatrix cMinusEqn
    (
    - fvm::laplacian(DMinus_nonD, cMinus)
    - fvm::div(DMinusphiE, cMinus)

    //- fvm::div(DMinus_nonD*zMinus*phiE, cMinus) 
    // - fvm::div(mobilityMinus * phiE, cMinus)
    //- fvc::laplacian(mobilityMinus * cMinus, psiE)
    );

    if(solveTransient)
    {
        cPlusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cPlus);
        cMinusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cMinus);
    }
    PNPEqn.insertEquation(1, cPlusEqn);
    PNPEqn.insertEquation(2, cMinusEqn);

    dimensionedScalar zPlusBy2sqrDebL_nonD("zPlusBy2sqrDebL_nonD", zPlus*oneOverSquareDebL_nonD);
    dimensionedScalar zMinusBy2sqrDebL_nonD("zMinusBy2sqrDebL_nonD", zMinus*oneOverSquareDebL_nonD);

    fvScalarMatrix Sp_fvScalarMatrix_cPlus(fvm::Sp(zPlusBy2sqrDebL_nonD, cPlus)); // (fvm::Sp((F*zPlus)/(eps0*epsr*NA), cPlus));
    fvScalarMatrix Sp_fvScalarMatrix_cMinus(fvm::Sp(zMinusBy2sqrDebL_nonD, cMinus)); // (fvm::Sp((F*zMinus)/(eps0*epsr*NA), cMinus));
    PNPEqn.insertEquationCoupling(0, 1, Sp_fvScalarMatrix_cPlus);
    PNPEqn.insertEquationCoupling(0, 2, Sp_fvScalarMatrix_cMinus);

    maxResidual = cmptMax(PNPEqn.solve().initialResidual()); // residual = b - A*x_n
    
    // Retrieve solution
    PNPEqn.retrieveSolution(0, psiE.internalField());
    PNPEqn.retrieveSolution(1, cPlus.internalField());
    PNPEqn.retrieveSolution(2, cMinus.internalField());

    psiE.correctBoundaryConditions();
    cPlus.correctBoundaryConditions();
    cMinus.correctBoundaryConditions();

    #include "boundC.H"
    #include "convergenceCheck.H"

}
