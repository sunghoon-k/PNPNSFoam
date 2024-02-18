netCharge = F*(zPlus*cPlus + zMinus*cMinus);
fvScalarMatrix psiEEqn
(
    fvm::laplacian(psiE)
    ==
    - netCharge/(eps0*epsr)
);
if (currentIter != nOuterIter -1)
{
    psiEEqn.relax();
}

// psiEEqn.solve();
psiEResidual = psiEEqn.solve().initialResidual(); // residual = b - A*x_n
psiE.correctBoundaryConditions();
