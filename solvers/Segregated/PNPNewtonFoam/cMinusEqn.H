{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  surfaceScalarField jMinus = -DMinus * (fvc::snGrad(cMinus) + (zMinus/psiE0) * fvc::interpolate(cMinus) * fvc::snGrad(psiE)) * mesh.magSf();
  fvScalarMatrix cMinusEqn
  (
  //    fvm::ddt(cPlus)
    - fvm::laplacian(DMinus, cMinus)
    //- fvm::div(DPlus_nonD * phiE, cMinus)
    //- fvm::div(mobilityPlus * phiE, cMinus)
    - fvc::laplacian((DMinus/psiE0) * cMinus, psiE)
    ==
    - fvc::surfaceIntegrate(jMinus)
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