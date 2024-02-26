{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /* Poisson Equation */

    // Add off-diagonal coupling terms
    fvScalarMatrix Sp_fvScalarMatrix_cPlus(fvm::Sp((e*zPlus)/(eps0*epsr), cPlus));
    fvScalarMatrix Sp_fvScalarMatrix_cMinus(fvm::Sp((e*zMinus)/(eps0*epsr), cMinus));

    PNPEqn.insertEquationCoupling(0, 1, Sp_fvScalarMatrix_cPlus); // 4: psiE <- 5: cPlus
    PNPEqn.insertEquationCoupling(0, 2, Sp_fvScalarMatrix_cMinus); // 4: psiE <- 6: cMinus
}
