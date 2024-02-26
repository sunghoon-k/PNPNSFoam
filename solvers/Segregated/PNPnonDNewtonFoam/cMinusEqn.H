{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  surfaceScalarField jMinus = -DMinus_nonD * (fvc::snGrad(cMinus) + zMinus * fvc::interpolate(cMinus) * fvc::snGrad(psiE)) * mesh.magSf();
  fvScalarMatrix cMinusEqn
  (
  //    fvm::ddt(cPlus)
    - fvm::laplacian(DMinus_nonD, cMinus)
    //- fvm::div(DPlus_nonD * phiE, cMinus)
    //- fvm::div(mobilityPlus * phiE, cMinus)
    - fvc::laplacian(DMinus_nonD * cMinus, psiE)
    ==
    - fvc::surfaceIntegrate(jMinus)
  );

  if(solveTransient)
  {
    cMinusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cMinus);
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