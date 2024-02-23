{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /* NS Equation */

    // Calculate grad p coupling matrix. Needs to be here if one uses
    // gradient schemes with limiters.  VV, 9/June/2014
    BlockLduSystem<vector, vector> pInU(fvm::grad(p));
    NSEqn.insertBlockCoupling(0, 3, pInU, true); // 0: U <- 3: p

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /* Continuity Equation */

    // Calculate div U coupling.  Could be calculated only once since
    // it is only geometry dependent.  VV, 9/June/2014
    BlockLduSystem<vector, scalar> UInp(fvm::UDiv(U));
    NSEqn.insertBlockCoupling(3, 0, UInp, false); // 3: p <- 0: U

}
