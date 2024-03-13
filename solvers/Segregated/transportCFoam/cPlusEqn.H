{
  //volScalarField mobilityPluscPlus(mobilityPlus * cPlus);
  surfaceScalarField mobilityPluscPlus(mobilityPlus * fvc::interpolate(cPlus));
  fvScalarMatrix cPlusEqn
  (
  //    fvm::ddt(cPlus)
    - fvm::laplacian(DPlus, cPlus)
    //- fvm::div(mobilityPlus * phiE, cPlus)
    // - fvc::laplacian(mobilityPlus * cPlus, psiE)
  );

  if(solveTransient)
  {
    cPlusEqn += fvm::ddt(cPlus);
  }

  if (currentIter == nOuterIter -1)
  {
    cPlusEqn.relax(1);
  }
  else
  {
    cPlusEqn.relax();
  }

  cPlusResidual = cPlusEqn.solve().initialResidual();
  cPlus.correctBoundaryConditions();

}