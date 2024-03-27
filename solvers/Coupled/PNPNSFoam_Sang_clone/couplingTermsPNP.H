// Prepare block system
fvBlockMatrix<vector3> PNPEqn(PNP);

fvScalarMatrix electroMigrationLaplacianPlusE
(
    fvm::laplacian(DPluszPluscPlusf, psiE)
);

fvScalarMatrix electroMigrationLaplacianMinusE
(
    fvm::laplacian(DMinuszMinuscMinusf, psiE)
);

PNPEqn.insertEquation(0, cPlusEqn);
PNPEqn.insertEquation(1, cMinusEqn);
PNPEqn.insertEquation(2, psiEEqn);

scalarField psiEcouplingPlus(mesh.nCells(),zPlus.value()*oneOver2SquareDebL_nonD);
scalarField psiEcouplingMinus(mesh.nCells(),zMinus.value()*oneOver2SquareDebL_nonD);

PNPEqn.insertEquationCoupling(0, 2, electroMigrationLaplacianPlusE); // 0: psiE -> 1: cPlus
PNPEqn.insertEquationCoupling(1, 2, electroMigrationLaplacianMinusE); // 0: psiE -> 2: cMinus
PNPEqn.insertEquationCoupling(2, 0, psiEcouplingPlus); // 0: psiE <- 1: cPlus
PNPEqn.insertEquationCoupling(2, 1, psiEcouplingMinus); // 0: psiE <- 2: cMinus
