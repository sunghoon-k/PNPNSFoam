{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  fvScalarMatrix cPlusEqn
  (
  //    fvm::ddt(cPlus)
    - fvm::laplacian(DPlus * NA, cPlus)
    - fvm::div(mobilityPlus * phiE * NA, cPlus)
    //- fvm::div(mobilityPlus * phiE, cPlus)
    //- fvc::laplacian(mobilityPlus * cPlus, psiE)
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

  //cPlusEqn *= (NA).value();

  cPlusResidual = cPlusEqn.solve().initialResidual();
  cPlus.correctBoundaryConditions();

}