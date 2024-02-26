{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /* Poisson Equation */

    // Add off-diagonal coupling terms
//    scalarField Phicoupling1(mesh.nCells(),(zPlus/(2*sqr(DebL_nonD))).value());
//    scalarField Phicoupling2(mesh.nCells(),(zMinus/(2*sqr(DebL_nonD))).value());

//    PNPEqn.insertEquationCoupling(0, 1, Phicoupling1);
//    PNPEqn.insertEquationCoupling(0, 2, Phicoupling2);
        
//    fvScalarMatrix Sp_fvScalarMatrix_cPlus(fvm::Sp(zPlus/(2*sqr(DebL_nonD)), cPlus)); // (fvm::Sp((F*zPlus)/(eps0*epsr*NA), cPlus));
//    fvScalarMatrix Sp_fvScalarMatrix_cMinus(fvm::Sp(zMinus/(2*sqr(DebL_nonD)), cMinus)); // (fvm::Sp((F*zMinus)/(eps0*epsr*NA), cMinus));
    dimensionedScalar zPlusBy2sqrDebL_nonD("zPlusBy2sqrDebL_nonD", zPlus/2*sqr(DebL_nonD));
    dimensionedScalar zMinusBy2sqrDebL_nonD("zMinusBy2sqrDebL_nonD", zMinus/2*sqr(DebL_nonD));

    fvScalarMatrix Sp_fvScalarMatrix_cPlus(fvm::Sp(zPlusBy2sqrDebL_nonD, cPlus)); // (fvm::Sp((F*zPlus)/(eps0*epsr*NA), cPlus));
    fvScalarMatrix Sp_fvScalarMatrix_cMinus(fvm::Sp(zMinusBy2sqrDebL_nonD, cMinus)); // (fvm::Sp((F*zMinus)/(eps0*epsr*NA), cMinus));
    PNPEqn.insertEquationCoupling(0, 1, Sp_fvScalarMatrix_cPlus);
    PNPEqn.insertEquationCoupling(0, 2, Sp_fvScalarMatrix_cMinus);


/*
    //- fvc::laplacian(mobilityPlus * cPlus, psiE)
    phicPlus = fvc::snGrad(cPlus) * mesh.magSf();
    phicMinus = fvc::snGrad(cMinus) * mesh.magSf();
    fvScalarMatrix psiEincPlus(- fvm::laplacian(DPlus_nonD * phicPlus, psiE));
    fvScalarMatrix psiEincMinus(- fvm::laplacian(DMinus_nonD * phicMinus, psiE));

    PNPEqn.insertEquationCoupling(1, 0, psiEincPlus); // 1: cPlus <- 0: psiE
    PNPEqn.insertEquationCoupling(2, 0, psiEincMinus); // 2: cMinus <- 0: psiE
*/

}
