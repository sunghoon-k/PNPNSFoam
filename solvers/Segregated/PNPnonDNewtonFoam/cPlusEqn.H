{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  surfaceScalarField J1 = -DPlus_nonD * (fvc::snGrad(cPlus) + zPlus * fvc::interpolate(cPlus) * fvc::snGrad(psiE)) * mesh.magSf();
  fvScalarMatrix cPlusEqn
  (
  //    fvm::ddt(cPlus)
    - fvm::laplacian(DPlus_nonD, cPlus)
    //- fvm::div(DPlus_nonD * phiE, cPlus)
    //- fvm::div(mobilityPlus * phiE, cPlus)
    - fvc::laplacian(DPlus_nonD * cPlus, psiE)
    ==
    - fvc::surfaceIntegrate(J1)
  );

  if(solveTransient)
  {
    cPlusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cPlus);
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