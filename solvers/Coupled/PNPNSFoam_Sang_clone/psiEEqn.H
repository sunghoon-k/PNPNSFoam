fvScalarMatrix psiEEqn
(
    fvm::laplacian(psiE)
);
psiEEqn.relax();

fvScalarMatrix psiIEqn
(
    fvm::laplacian(psiI)
);

psiIEqn.relax();
