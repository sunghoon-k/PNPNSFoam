fvScalarMatrix psiEEqn
(
    fvm::laplacian(psiE)
  //+ (zPlus*cPlus+zMinus*cMinus)/sqr(l0_one) 
  //+ F3
//  - fvc::laplacian(psiE)
//  - (zPlus*cPlus+zMinus*cMinus)/sqr(l0_one)
);
psiEEqn.relax();
psiEEqn.setReference(0, 0);

/*
fvScalarMatrix psiIEqn
(
    fvm::laplacian(dpsiI)
  + fvc::laplacian(psiI)
);

psiIEqn.relax();

*/
