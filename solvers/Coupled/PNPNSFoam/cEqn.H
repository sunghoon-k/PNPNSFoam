{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  fvScalarMatrix cPlusEqn
  (
      fvm::div(phi, cPlus)
    - fvm::laplacian(DPlus, cPlus)
    - fvm::div(mobilityPlus * phiE, cPlus)
    //- fvm::div(mobilityPlus * phiE, cPlus)
    //- fvc::laplacian(mobilityPlus * cPlus, psiE)
  );

  fvScalarMatrix cMinusEqn
  (
      fvm::div(phi, cMinus)
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
  PNPNSEqn.insertEquation(5, cPlusEqn);
  PNPNSEqn.insertEquation(6, cMinusEqn);

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