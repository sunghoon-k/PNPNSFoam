{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /* NS Equation */

    // Calculate grad p coupling matrix. Needs to be here if one uses
    // gradient schemes with limiters.  VV, 9/June/2014
    BlockLduSystem<vector, vector> pInU(fvm::grad(p));
    pInU *= UScaling; // BlockLduSystem 에는 scalar.value()를 곱해야됨.

    volScalarField rhoE( psiE * dimensionedScalar("norm", eps0.dimensions()/dimArea, 0.) );
    rhoE += (zPlus*cPlus*F);
    rhoE += (zMinus*cMinus*F);

    BlockLduSystem<vector, vector> psiEInU(fvm::grad(psiE));
    // rhoE/rho;
    psiEInU *= (rhoE.internalField()/rho.value()) * UScaling ;

    PNPNSEqn.insertBlockCoupling(0, 3, pInU, true); // 0: U <- 3: p
    PNPNSEqn.insertBlockCoupling(0, 4, psiEInU, true); // 0: U <- 4: psiE

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /* Continuity Equation */

    // Calculate div U coupling.  Could be calculated only once since
    // it is only geometry dependent.  VV, 9/June/2014
    BlockLduSystem<vector, scalar> UInp(fvm::UDiv(U));

    UInp *= pScaling;
    PNPNSEqn.insertBlockCoupling(3, 0, UInp, false); // 3: p <- 0: U

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /* Poisson Equation */

    // Add off-diagonal coupling terms
    fvScalarMatrix Sp_fvScalarMatrix_cPlus(fvm::Sp((F*zPlus)/(eps0*epsr), cPlus));
    fvScalarMatrix Sp_fvScalarMatrix_cMinus(fvm::Sp((F*zMinus)/(eps0*epsr), cMinus));

    // fvScalarMatrix SuSp_fvScalarMatrix_cPlus(fvm::SuSp(cPlus));
    // fvScalarMatrix SuSp_fvScalarMatrix_cMinus(fvm::SuSp(cMinus));

    Sp_fvScalarMatrix_cPlus *= psiEScaling;
    Sp_fvScalarMatrix_cMinus *= psiEScaling;

    PNPNSEqn.insertEquationCoupling(4, 5, Sp_fvScalarMatrix_cPlus); // 4: psiE <- 5: cPlus
    PNPNSEqn.insertEquationCoupling(4, 6, Sp_fvScalarMatrix_cMinus); // 4: psiE <- 6: cMinus
}
