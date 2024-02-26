{
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
    cPlusEqn += fvm::ddt(cPlus) + fvc::ddt(cPlus_old);
    cMinusEqn += fvm::ddt(cMinus) + fvc::ddt(cMinus_old);
  }

  PNPEqn.insertEquation(1, cPlusEqn);
  PNPEqn.insertEquation(2, cMinusEqn);
}