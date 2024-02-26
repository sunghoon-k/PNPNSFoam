{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  fvScalarMatrix cMinusEqn
  (
  //    fvm::ddt(cPlus)
    - fvm::laplacian(DMinus_nonD, cMinus)
    - fvm::div(DMinus_nonD * phiE, cMinus)

    //- fvm::div(mobilityMinus * phiE, cMinus)
    //- fvc::laplacian(mobilityMinus * cMinus, psiE)
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