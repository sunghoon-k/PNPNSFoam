{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  surfaceScalarField J1 = -DPlus * (fvc::snGrad(cPlus) + (zPlus/psiE0) * fvc::interpolate(cPlus) * fvc::snGrad(psiE)) * mesh.magSf();
  fvScalarMatrix cPlusEqn
  (
  //    fvm::ddt(cPlus)
    - fvm::laplacian(DPlus, cPlus)
    //- fvm::div(DPlus_nonD * phiE, cPlus)
    //- fvm::div(mobilityPlus * phiE, cPlus)
    - fvc::laplacian((DPlus/psiE0) * cPlus, psiE)
    ==
    - fvc::surfaceIntegrate(J1)
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