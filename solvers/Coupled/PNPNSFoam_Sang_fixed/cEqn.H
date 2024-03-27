fvScalarMatrix cPlusEqn
(
    fvm::laplacian(-DPlus_nonD, cPlus)
  + fvm::div(DPluszPlusphiE,    cPlus)
//  - fvc::laplacian(-DPlus_nonD, cPlus)
//  - fvc::div(DPluszPlusphiE, cPlus)
  //+ F1
);

fvScalarMatrix cMinusEqn
(
    fvm::laplacian(-DMinus_nonD,  cMinus)
  + fvm::div(DMinuszMinusphiE,    cMinus)
//  - fvc::laplacian(-DMinus_nonD, cMinus)
//  - fvc::div(DMinuszMinusphiE, cMinus)
  //+ F2
);

if(solveTransient)
{
  cPlusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cPlus);
  cMinusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cMinus);
}

cPlusEqn.relax();
cMinusEqn.relax();

cPlusEqn.setReference(0, 0);
cMinusEqn.setReference(0, 0);
