fvScalarMatrix electroMigrationLaplacianPlusI
(
    fvm::laplacian(DPluszPluscPlusf, psiI)
);

fvScalarMatrix electroMigrationLaplacianMinusI
(
    fvm::laplacian(DMinuszMinuscMinusf, psiI)
);

fvScalarMatrix electroMigrationLaplacianPlusE
(
    fvm::laplacian(DPluszPluscPlusf, psiE)
);

fvScalarMatrix electroMigrationLaplacianMinusE
(
    fvm::laplacian(DMinuszMinuscMinusf, psiE)
);

PNPNSsEqn.insertEquation(0, cPlusEqn);
PNPNSsEqn.insertEquation(1, cMinusEqn);
PNPNSsEqn.insertEquation(2, psiIEqn);
PNPNSsEqn.insertEquation(3, psiEEqn);
PNPNSsEqn.insertEquation(4, UEqn);
PNPNSsEqn.insertEquation(7, pEqn);

scalarField psiEcouplingPlus(mesh.nCells(),zPlus.value()*oneOver2SquareDebL_nonD);
scalarField psiEcouplingMinus(mesh.nCells(),zMinus.value()*oneOver2SquareDebL_nonD);

PNPNSsEqn.insertEquationCoupling(0, 2, electroMigrationLaplacianPlusI); // 0: psiE -> 1: cPlus
PNPNSsEqn.insertEquationCoupling(1, 2, electroMigrationLaplacianMinusI); // 0: psiE -> 2: cMinus
PNPNSsEqn.insertEquationCoupling(0, 3, electroMigrationLaplacianPlusE); // 0: psiE -> 1: cPlus
PNPNSsEqn.insertEquationCoupling(1, 3, electroMigrationLaplacianMinusE); // 0: psiE -> 2: cMinus
PNPNSsEqn.insertEquationCoupling(2, 0, psiEcouplingPlus); // 0: psiE <- 1: cPlus
PNPNSsEqn.insertEquationCoupling(2, 1, psiEcouplingMinus); // 0: psiE <- 2: cMinus

BlockLduSystem<vector, vector> pInU(fvm::grad(p));
PNPNSsEqn.insertBlockCoupling(4, 7, pInU, true); // 0: U <- 3: p

BlockLduSystem<vector, scalar> UInp(fvm::UDiv(U));
PNPNSsEqn.insertBlockCoupling(7, 4, UInp, false); // 3: p <- 0: U
