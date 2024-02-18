{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  fvScalarMatrix cPlusEqn
  (
    - fvm::laplacian(DPlus, cPlus)
    - fvm::div(mobilityPlus * phiE, cPlus)
    //- fvm::div(mobilityPlus * phiE, cPlus)
    //- fvc::laplacian(mobilityPlus * cPlus, psiE)
  );

  fvScalarMatrix cMinusEqn
  (
    - fvm::laplacian(DMinus, cMinus)
    - fvm::div(mobilityMinus * phiE, cMinus)
    //- fvm::div(mobilityMinus * phiE, cMinus)
    //- fvc::laplacian(mobilityMinus * cMinus, psiE)
  );

  if(solveTransient)
  {
    cPlusEqn += fvm::ddt(cPlus);
    cMinusEqn += fvm::ddt(cMinus);
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