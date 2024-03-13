netCharge = e*(zPlus*cPlus + zMinus*cMinus);
// Info << netCharge.dimensions() << endl;
volScalarField netChargeOverEpsilon(netCharge/(eps0*epsr));

fvScalarMatrix psiEEqn
(
    fvm::laplacian(psiE)
    ==
    - netChargeOverEpsilon
);
if (currentIter != nOuterIter -1)
{
    psiEEqn.relax();
}

// psiEEqn.solve();
//psiEEqn.source() *= (1/psiE0).value();
psiEResidual = psiEEqn.solve().initialResidual(); // residual = b - A*x_n

//psiE.internalField() *= psiE0.value();
psiE.correctBoundaryConditions();
