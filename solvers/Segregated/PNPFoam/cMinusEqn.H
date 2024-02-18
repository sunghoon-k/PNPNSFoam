{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  fvScalarMatrix cMinusEqn
  (
  //    fvm::ddt(cPlus)
    - fvm::laplacian(DMinus*NA, cMinus)
    - fvm::div(mobilityMinus*NA * phiE, cMinus)

    //- fvm::div(mobilityMinus * phiE, cMinus)
    //- fvc::laplacian(mobilityMinus * cMinus, psiE)
  );

  if(solveTransient)
  {
    cMinusEqn += fvm::ddt(cMinus);
  }

  if (currentIter == nOuterIter -1)
  {
    cMinusEqn.relax(1);
  }
  else
  {
    cMinusEqn.relax();
  }

  //cMinusEqn *= (NA).value();
  cMinusResidual = cMinusEqn.solve().initialResidual();

  //Info << "calculate boundary condition" << nl << endl;
  cMinus.correctBoundaryConditions();
}