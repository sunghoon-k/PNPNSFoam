fvScalarMatrix cPlusEqn
(
    fvm::laplacian(-DPlus_nonD, cPlus)
  + fvm::div(DPluszPlusphiE, cPlus)
);

fvScalarMatrix cMinusEqn
(
    fvm::laplacian(-DMinus_nonD, cMinus)
  + fvm::div(DMinuszMinusphiE, cMinus)
);

if(solveTransient)
{
  cPlusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cPlus);
  cMinusEqn += (t0_one/sqr(l0_one)) * fvm::ddt(cMinus);
}

cPlusEqn.relax();
cMinusEqn.relax();
