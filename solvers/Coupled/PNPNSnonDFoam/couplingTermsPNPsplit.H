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

PNPsEqn.insertEquation(0, cPlusEqn);
PNPsEqn.insertEquation(1, cMinusEqn);
PNPsEqn.insertEquation(2, psiIEqn);
PNPsEqn.insertEquation(3, psiEEqn);

scalarField psiEcouplingPlus(mesh.nCells(),zPlus.value()*oneOver2SquareDebL_nonD);
scalarField psiEcouplingMinus(mesh.nCells(),zMinus.value()*oneOver2SquareDebL_nonD);

PNPsEqn.insertEquationCoupling(0, 2, electroMigrationLaplacianPlusI); // 0: psiE -> 1: cPlus
PNPsEqn.insertEquationCoupling(1, 2, electroMigrationLaplacianMinusI); // 0: psiE -> 2: cMinus
PNPsEqn.insertEquationCoupling(0, 3, electroMigrationLaplacianPlusE); // 0: psiE -> 1: cPlus
PNPsEqn.insertEquationCoupling(1, 3, electroMigrationLaplacianMinusE); // 0: psiE -> 2: cMinus
PNPsEqn.insertEquationCoupling(2, 0, psiEcouplingPlus); // 0: psiE <- 1: cPlus
PNPsEqn.insertEquationCoupling(2, 1, psiEcouplingMinus); // 0: psiE <- 2: cMinus