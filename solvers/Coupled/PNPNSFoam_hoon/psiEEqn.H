fvScalarMatrix psiEEqn
(
    fvm::laplacian(psiE)
);
psiEEqn.relax();
psiEEqn.setReference(0, 0);

/*
fvScalarMatrix psiIEqn
(
    fvm::laplacian(dpsiI)
  + fvc::laplacian(psiI)
);

psiIEqn.relax();

*/
