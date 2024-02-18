{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  fvScalarMatrix cMinusEqn
  (
  //    fvm::ddt(cPlus)
    - fvm::laplacian(DMinus, cMinus)
    - fvm::div(mobilityMinus * phiE, cMinus)
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

  cMinusResidual = cMinusEqn.solve().initialResidual();
  cMinus.correctBoundaryConditions();
}