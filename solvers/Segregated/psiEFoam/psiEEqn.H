
fvScalarMatrix psiEEqn
(
    fvm::laplacian(psiE)
);

psiEEqn.relax();
psiE *= psiEScaling;
psiEEqn.solve();
maxResidual = cmptMax(PNPNSEqn.solve().finalResidual()); // residual = b - A*x_n
