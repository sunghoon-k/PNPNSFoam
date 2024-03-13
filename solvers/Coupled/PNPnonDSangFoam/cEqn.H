fvScalarMatrix cPlusEqn
(
  fvm::laplacian(-DPlus_nonD, cPlus)
+ fvm::div(DPluszPlusphiE, cPlus)
//- fvm::div(mobilityPlus * phiE, cPlus)
//- fvc::laplacian(DPluszPluscPlusf, psiE)
//+ fvc::div(mobilityPlus * phiE, cPlus)
//==
//- FcPlus
);
cPlusEqn.relax();

fvScalarMatrix cMinusEqn
(
  fvm::laplacian(-DMinus_nonD, cMinus)
+ fvm::div(DMinuszMinusphiE, cMinus)
//- fvm::div(mobilityMinus * phiE, cMinus)
//- fvc::laplacian(DMinuszMinuscMinusf, psiE)
//+ fvc::div(mobilityPlus * phiE, cPlus)
//==
//- FcMinus
);
cMinusEqn.relax();

if(solveTransient)
{
cPlusEqn += fvm::ddt(cPlus);
cMinusEqn += fvm::ddt(cMinus);
}