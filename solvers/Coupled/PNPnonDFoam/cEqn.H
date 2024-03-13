
fvScalarMatrix cPlusEqn
(
    fvm::laplacian(-DPlus_nonD, cPlus)
  + fvm::div(DPluszPlusphiE, cPlus)
  //+ fvm::div(DPluszPlusinterphiE, cPlus)
  //+ fvc::laplacian(-DPlus_nonD, cPlus)
  //+ fvc::div(DPluszPlusphiE, cPlus)
  //- fvc::surfaceIntegrate(jPlus) 
  //- fvm::div(DPlus_nonD*zPlus*phiE, cPlus) // 
  //- fvm::div(mobilityPlus * phiE, cPlus)
  //- fvc::laplacian(phiDcPlus, psiE)
);

fvScalarMatrix cMinusEqn
(
    fvm::laplacian(-DMinus_nonD, cMinus)
  + fvm::div(DMinuszMinusphiE, cMinus)
  //+ fvm::div(DMinuszMinusinterphiE, cMinus)
  //+ fvc::laplacian(-DMinus_nonD, cMinus)
  //+ fvc::div(DMinuszMinusphiE, cMinus)
  //- fvc::surfaceIntegrate(jMinus) 

  //- fvm::div(DMinus_nonD*zMinus*phiE, cMinus) 
  // - fvm::div(mobilityMinus * phiE, cMinus)
  //- fvc::laplacian(phiDcMinus, psiE)
);

if(solveTransient)
{
  cPlusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cPlus);
  cMinusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cMinus);
}

cPlusEqn.relax();
cMinusEqn.relax();


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
