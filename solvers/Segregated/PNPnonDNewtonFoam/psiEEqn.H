netCharge = (zPlus*cPlus + zMinus*cMinus);
// Info << netCharge.dimensions() << endl;
surfaceScalarField jpsiE = fvc::snGrad(psiE);

fvScalarMatrix psiEEqn
(
    fvm::laplacian(psiE)
  + netCharge/sqr(DebL_nonD * l0_one)
  ==
  - fvc::surfaceIntegrate(jpsiE)
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
