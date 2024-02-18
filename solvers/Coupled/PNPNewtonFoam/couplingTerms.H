{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /* Poisson Equation */

    // Add off-diagonal coupling terms
    fvScalarMatrix Sp_fvScalarMatrix_cPlus(fvm::Sp((F*zPlus)/(eps0*epsr), cPlus));
    fvScalarMatrix Sp_fvScalarMatrix_cMinus(fvm::Sp((F*zMinus)/(eps0*epsr), cMinus));

    // fvScalarMatrix SuSp_fvScalarMatrix_cPlus(fvm::SuSp(cPlus));
    // fvScalarMatrix SuSp_fvScalarMatrix_cMinus(fvm::SuSp(cMinus));

    Sp_fvScalarMatrix_cPlus *= psiEScaling;
    Sp_fvScalarMatrix_cMinus *= psiEScaling;

    PNPEqn.insertEquationCoupling(0, 1, Sp_fvScalarMatrix_cPlus); // 4: psiE <- 5: cPlus
    PNPEqn.insertEquationCoupling(0, 2, Sp_fvScalarMatrix_cMinus); // 4: psiE <- 6: cMinus
}
