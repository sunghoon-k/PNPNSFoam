{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  fvScalarMatrix cPlusEqn
  (
    - fvm::laplacian(DPlus_nonD, cPlus)
    //- fvm::div(DPlus_nonD*zPlus*phiE, cPlus) // - fvm::div(mobilityPlus * phiE, cPlus)
    //- fvc::laplacian(mobilityPlus * cPlus, psiE)
  );

  fvScalarMatrix cMinusEqn
  (
    - fvm::laplacian(DMinus_nonD, cMinus)
    //- fvm::div(DMinus_nonD*zMinus*phiE, cMinus) // - fvm::div(mobilityMinus * phiE, cMinus)
    //- fvc::laplacian(mobilityMinus * cMinus, psiE)
  );

  if(solveTransient)
  {
    cPlusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cPlus);
    cMinusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cMinus);
  }
  PNPEqn.insertEquation(1, cPlusEqn);
  PNPEqn.insertEquation(2, cMinusEqn);

/*
  if (currentIter == nOuterIter -1)
  {
    cPlusEqn.relax(1);
  }
  else
  {
    cPlusEqn.relax();
  }
*/

}