{
  // surfaceScalarField jcPlus = -DPlus * (fvc::snGrad(cPlus_old) + (zPlus/psiE0) * fvc::interpolate(cPlus_old) * fvc::snGrad(psiE_old));
  fvScalarMatrix cPlusEqn
  (
  //    fvm::ddt(cPlus)
    - fvm::laplacian(DPlus, cPlus)
    - fvm::div((DPlus/psiE0) * phiE_old, cPlus)
    //- fvm::div(mobilityPlus * phiE, cPlus)
    // - fvc::laplacian((DPlus/psiE0) * cPlus, psiE)
    ==
    - ( - fvc::laplacian(DPlus, cPlus_old) - fvc::div((DPlus/psiE0) * phiE_old, cPlus_old) )
  );

  if(solveTransient)
  {
    cPlusEqn += fvm::ddt(cPlus) + fvc::ddt(cPlus_old);
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
  cPlus += cPlus_old;
  cPlus.correctBoundaryConditions();
  cPlus_old = cPlus;
}