{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /* NS Equation */

    // Calculate grad p coupling matrix. Needs to be here if one uses
    // gradient schemes with limiters.  VV, 9/June/2014
    BlockLduSystem<vector, vector> pInU(fvm::grad(p));

    BlockLduSystem<vector, vector> psiEInU(fvm::grad(psiE));
    // rhoE/rho;
    volScalarField netChargeoneOverSquareDebL_nonD(netCharge * oneOverSquareDebL_nonD);
    psiEInU *= netChargeoneOverSquareDebL_nonD.internalField();

    PNPNSEqn.insertBlockCoupling(0, 3, pInU, true); // 0: U <- 3: p
    PNPNSEqn.insertBlockCoupling(0, 4, psiEInU, true); // 0: U <- 4: psiE

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /* Continuity Equation */

    // Calculate div U coupling.  Could be calculated only once since
    // it is only geometry dependent.  VV, 9/June/2014
    BlockLduSystem<vector, scalar> UInp(fvm::UDiv(U));

    PNPNSEqn.insertBlockCoupling(3, 0, UInp, false); // 3: p <- 0: U

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /* Poisson Equation */

    dimensionedScalar zPlusBy2sqrDebL_nonD("zPlusBy2sqrDebL_nonD", zPlus*oneOverSquareDebL_nonD);
    dimensionedScalar zMinusBy2sqrDebL_nonD("zMinusBy2sqrDebL_nonD", zMinus*oneOverSquareDebL_nonD);

    //oneField oneFiedl_temp();
    
    // fvm::SuSp(zPlusBy2sqrDebL_nonD, cPlus);
    fvScalarMatrix Sp_fvScalarMatrix_cPlus(fvm::Sp(zPlusBy2sqrDebL_nonD, cPlus)); // (fvm::Sp((F*zPlus)/(eps0*epsr*NA), cPlus));
    fvScalarMatrix Sp_fvScalarMatrix_cMinus(fvm::Sp(zMinusBy2sqrDebL_nonD, cMinus)); // (fvm::Sp((F*zMinus)/(eps0*epsr*NA), cMinus));
    PNPNSEqn.insertEquationCoupling(4, 5, Sp_fvScalarMatrix_cPlus);
    PNPNSEqn.insertEquationCoupling(4, 6, Sp_fvScalarMatrix_cMinus);
}
