//netCharge = e*(zPlus*cPlus + zMinus*cMinus);
//Info << netCharge.dimensions() << endl;
/*

*/
fvScalarMatrix psiEEqn
(
      fvm::laplacian(psiE)
    //==
    //  - FpsiE
);
psiEEqn.relax();

/*
if (currentIter != nOuterIter -1)
{
    psiEEqn.relax();
}
*/