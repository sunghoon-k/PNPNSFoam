{
  // surfaceScalarField jcMinus = -DMinus * (fvc::snGrad(cMinus_old) + (zMinus/psiE0) * fvc::interpolate(cMinus_old) * fvc::snGrad(psiE_old));
  fvScalarMatrix cMinusEqn
  (
  //    fvm::ddt(cMinus)
    - fvm::laplacian(DMinus, cMinus)
    - fvm::div((DMinus/psiE0) * phiE_old, cMinus)
    //- fvm::div(mobilityMinus * phiE, cMinus)
    // - fvc::laplacian((DMinus/psiE0) * cMinus, psiE)
    ==
    - ( - fvc::laplacian(DMinus, cMinus_old) - fvc::div((DMinus/psiE0) * phiE_old, cMinus_old) )
  );

  if(solveTransient)
  {
    cMinusEqn += fvm::ddt(cMinus) + fvc::ddt(cMinus_old);
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
  cMinus += cMinus_old;
  cMinus.correctBoundaryConditions();
  cMinus_old = cMinus;

}