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

/*
fvScalarMatrix electroMigrationLaplacianPlus2
(
    fvm::laplacian(DPluszPluscPlusf, psiI)
);
electroMigrationLaplacianPlus.relax();


fvScalarMatrix electroMigrationLaplacianMinus2
(
    fvm::laplacian(DMinuszMinuscMinusf, psiI)
);
electroMigrationLaplacianMinus.relax();

*/


PNPEqn.insertEquation(0, cPlusEqn);
PNPEqn.insertEquation(1, cMinusEqn);
PNPEqn.insertEquation(2, psiEEqn);
//PNPEqn.insertEquation(3, psiEEqn);

/*
fvScalarMatrix psiEcouplingPlus
(
    fvm::Sp(zPlus*sqrDebL_nonD/2, cPlus)
);
psiEcouplingPlus.relax();

fvScalarMatrix psiEcouplingMinus
(
    fvm::Sp(zPlus*sqrDebL_nonD/2, cMinus)
);
psiEcouplingMinus.relax();

*/
scalarField psiEcouplingPlus(mesh.nCells(),zPlus.value()*sqrDebL_nonD/2);
scalarField psiEcouplingMinus(mesh.nCells(),zMinus.value()*sqrDebL_nonD/2);

//PNPEqn.insertEquationCoupling(0, 3, electroMigrationLaplacianPlus); // 0: psiE -> 1: cPlus
//PNPEqn.insertEquationCoupling(1, 3, electroMigrationLaplacianMinus); // 0: psiE -> 2: cMinus
PNPEqn.insertEquationCoupling(0, 2, electroMigrationLaplacianPlus); // 0: psiE -> 1: cPlus
PNPEqn.insertEquationCoupling(1, 2, electroMigrationLaplacianMinus); // 0: psiE -> 2: cMinus
PNPEqn.insertEquationCoupling(2, 0, psiEcouplingPlus); // 0: psiE <- 1: cPlus
PNPEqn.insertEquationCoupling(2, 1, psiEcouplingMinus); // 0: psiE <- 2: cMinus

//PNPEqn.insertEquationCoupling(0, 2, electroMigrationLaplacianPlus2); // 0: psiE -> 1: cPlus
//PNPEqn.insertEquationCoupling(1, 2, electroMigrationLaplacianMinus2); // 0: psiE -> 2: cMinus
