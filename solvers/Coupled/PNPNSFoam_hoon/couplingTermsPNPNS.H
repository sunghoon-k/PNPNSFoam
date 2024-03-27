// Prepare block system
fvBlockMatrix<vector7> PNPNSEqn(PNPNS);

fvScalarMatrix electroMigrationLaplacianPlus
(
    fvm::laplacian(DPluszPluscPlusf, psiE)
);
electroMigrationLaplacianPlus.relax();


fvScalarMatrix electroMigrationLaplacianMinus
(
    fvm::laplacian(DMinuszMinuscMinusf, psiE)
);
electroMigrationLaplacianMinus.relax();

PNPNSEqn.insertEquation(0, cPlusEqn);
PNPNSEqn.insertEquation(1, cMinusEqn);
PNPNSEqn.insertEquation(2, psiEEqn);
PNPNSEqn.insertEquation(3, UEqn);
PNPNSEqn.insertEquation(6, pEqn);

scalarField psiEcouplingPlus(mesh.nCells(),zPlus.value()*oneOver2SquareDebL_nonD);
scalarField psiEcouplingMinus(mesh.nCells(),zMinus.value()*oneOver2SquareDebL_nonD);

PNPNSEqn.insertEquationCoupling(0, 2, electroMigrationLaplacianPlus); // 0: psiE -> 1: cPlus
PNPNSEqn.insertEquationCoupling(1, 2, electroMigrationLaplacianMinus); // 0: psiE -> 2: cMinus
PNPNSEqn.insertEquationCoupling(2, 0, psiEcouplingPlus); // 0: psiE <- 1: cPlus
PNPNSEqn.insertEquationCoupling(2, 1, psiEcouplingMinus); // 0: psiE <- 2: cMinus

BlockLduSystem<vector, vector> pInU(fvm::grad(p));
PNPNSEqn.insertBlockCoupling(3, 6, pInU, true); // 3: U <- 6: p

BlockLduSystem<vector, scalar> UInp(fvm::UDiv(U));
PNPNSEqn.insertBlockCoupling(6, 3, UInp, false); // 6: p <- 3: U


netCharge = zPlus*cPlus + zMinus*cMinus;

bodyForce = - netCharge * fvc::grad(psiE) * oneOver2SquareDebL_nonD;

BlockLduSystem<vector, vector> psiEInU(fvm::grad(psiE));
psiEInU *= (netCharge.internalField() * oneOver2SquareDebL_nonD);
PNPNSEqn.insertBlockCoupling(3, 2, psiEInU, true); // 0: U <- 2: psiE
