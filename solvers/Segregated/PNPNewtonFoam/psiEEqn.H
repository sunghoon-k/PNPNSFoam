netCharge = e*(zPlus*cPlus + zMinus*cMinus);
// Info << netCharge.dimensions() << endl;
// surfaceScalarField jpsiE = fvc::snGrad(psiE) * mesh.Sf();

fvScalarMatrix psiEEqn
(
    fvm::laplacian(psiE)
  + netCharge/(eps0*epsr)
  ==
  - fvc::laplacian(psiE_old)
);
if (currentIter != nOuterIter -1)
{
    psiEEqn.relax();
}

// psiEEqn.solve();
//psiEEqn.source() *= (1/psiE0).value();
psiEResidual = psiEEqn.solve().initialResidual(); // residual = b - A*x_n
psiE += psiE_old;

//psiE.internalField() *= psiE0.value();
psiE.correctBoundaryConditions();
psiE_old = psiE;
