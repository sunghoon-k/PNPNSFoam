{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /* NS Equation */

    // Calculate grad p coupling matrix. Needs to be here if one uses
    // gradient schemes with limiters.  VV, 9/June/2014
    BlockLduSystem<vector, vector> pInU(fvm::grad(p));

    BlockLduSystem<vector, vector> psiEInU(fvm::grad(psiE));
    // rhoE/rho;
    psiEInU *= (netCharge.internalField()/rho.value());

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

    // Add off-diagonal coupling terms
    fvScalarMatrix Sp_fvScalarMatrix_cPlus(fvm::SuSp((e*zPlus)/(eps0*epsr), cPlus));
    fvScalarMatrix Sp_fvScalarMatrix_cMinus(fvm::SuSp((e*zMinus)/(eps0*epsr), cMinus));

    PNPNSEqn.insertEquationCoupling(4, 5, Sp_fvScalarMatrix_cPlus); // 4: psiE <- 5: cPlus
    PNPNSEqn.insertEquationCoupling(4, 6, Sp_fvScalarMatrix_cMinus); // 4: psiE <- 6: cMinus
}
