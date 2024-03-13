{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  DPlusphiE = DPlus_nonD*zPlus*phiE;
  DMinusphiE = DMinus_nonD*zMinus*phiE;

  fvScalarMatrix cPlusEqn
  (
      Pe * fvm::div(phi, cPlus)
    - fvm::laplacian(DPlus_nonD, cPlus)
    - fvm::div(DPlusphiE, cPlus)
    //- fvm::div(DPlus_nonD*zPlus*phiE, cPlus) // 
    //- fvm::div(mobilityPlus * phiE, cPlus)
    //- fvc::laplacian(mobilityPlus * cPlus, psiE)
  );

  fvScalarMatrix cMinusEqn
  (
    Pe * fvm::div(phi, cMinus)
    - fvm::laplacian(DMinus_nonD, cMinus)
    - fvm::div(DMinusphiE, cMinus)

    //- fvm::div(DMinus_nonD*zMinus*phiE, cMinus) 
    // - fvm::div(mobilityMinus * phiE, cMinus)
    //- fvc::laplacian(mobilityMinus * cMinus, psiE)
  );

  if(solveTransient)
  {
    cPlusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cPlus);
    cMinusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cMinus);
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