p.storePrevIter();

#include "UEqn.H"
#include "pEqn.H"

BlockLduSystem<vector, vector> pInU(fvm::grad(p));
PNPNSEqn.insertBlockCoupling(3, 6, pInU, true); // 0: U <- 3: p

BlockLduSystem<vector, scalar> UInp(fvm::UDiv(U));
PNPNSEqn.insertBlockCoupling(6, 3, UInp, false); // 3: p <- 0: U


