PNPEqn.insertEquation(0, cPlusEqn);
PNPEqn.insertEquation(1, cMinusEqn);
PNPEqn.insertEquation(2, psiEEqn);

fvScalarMatrix electroMigrationLaplacianPlus(fvm::laplacian(-DPluszPluscPlusf, psiE));
fvScalarMatrix electroMigrationLaplacianMinus(fvm::laplacian(-DMinuszMinuscMinusf, psiE));
PNPEqn.insertEquationCoupling(0, 2, electroMigrationLaplacianPlus); // 0: psiE -> 1: cPlus
PNPEqn.insertEquationCoupling(1, 2, electroMigrationLaplacianMinus); // 0: psiE -> 2: cMinus


scalarField psiEcouplingPlus(mesh.nCells(),(zPlus.value()/(2*sqr(DebL_nonD.value()))));
scalarField psiEcouplingMinus(mesh.nCells(),(zMinus.value()/(2*sqr(DebL_nonD.value()))));
PNPEqn.insertEquationCoupling(2, 0, psiEcouplingPlus); // 0: psiE <- 1: cPlus
PNPEqn.insertEquationCoupling(2, 1, psiEcouplingMinus); // 0: psiE <- 2: cMinus

