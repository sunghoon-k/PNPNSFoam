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

//volScalarField temporarycPlus = ddtc_one * fvc::ddt(cPlus);
//Info << temporarycPlus.dimensions() << endl;
//Info << fvc::ddt(cPlus).dimensions() << endl;
//(t0_one/sqr(l0_one)) * 

if(pseudoTransient or solveTransient)
{
  cPlusEqn += ddtc_one * fvm::ddt(cPlus);
  cMinusEqn += ddtc_one * fvm::ddt(cMinus);
}


cPlusEqn.relax();
cMinusEqn.relax();

cPlusEqn.setReference(0, 0);
cMinusEqn.setReference(0, 0);
