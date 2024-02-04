{
    // Calculate grad p coupling matrix. Needs to be here if one uses
    // gradient schemes with limiters.  VV, 9/June/2014
    BlockLduSystem<vector, vector> pInU(fvm::grad(p));

    volScalarField rhoE( phiE * dimensionedScalar("norm", epsilon0.dimensions()/dimArea, 0.) );
    rhoE += (zPlus*cPlus*F_const);
    rhoE += (zMinus*cMinus*F_const);

    BlockLduSystem<vector, vector> phiEInU(rhoE/rho*fvm::grad(phiE));

    // Calculate div U coupling.  Could be calculated only once since
    // it is only geometry dependent.  VV, 9/June/2014
    BlockLduSystem<vector, scalar> UInp(fvm::UDiv(U));

    // Last argument in insertBlockCoupling says if the column direction
    // should be incremented. This is needed for arbitrary positioning
    // of U and p in the system. This could be better. VV, 30/April/2014
    PNPNSEqn.insertBlockCoupling(0, 3, pInU, true);
    PNPNSEqn.insertBlockCoupling(0, 4, phiEInU, true);

    PNPNSEqn.insertBlockCoupling(3, 0, UInp, false);
}
