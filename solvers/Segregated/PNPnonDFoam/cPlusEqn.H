{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  fvScalarMatrix cPlusEqn
  (
  //    fvm::ddt(cPlus)
    - fvm::laplacian(DPlus_nonD, cPlus)
    - fvm::div(DPlus_nonD * phiE, cPlus)
    //- fvm::div(mobilityPlus * phiE, cPlus)
    //- fvc::laplacian(mobilityPlus * cPlus, psiE)
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