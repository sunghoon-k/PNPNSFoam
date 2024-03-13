{

//  volScalarField mobilityMinuscMinus(mobilityMinus * cMinus);
  surfaceScalarField mobilityMinuscMinus(mobilityMinus * fvc::interpolate(cMinus));

  fvScalarMatrix cMinusEqn
  (
  //    fvm::ddt(cPlus)
    - fvm::laplacian(DMinus, cMinus)
    //- fvm::div(mobilityMinus * phiE, cMinus)

    //- fvm::div(mobilityMinus * phiE, cMinus)
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

  cMinusResidual = cMinusEqn.solve().initialResidual();

  //Info << "calculate boundary condition" << nl << endl;
  cMinus.correctBoundaryConditions();
}