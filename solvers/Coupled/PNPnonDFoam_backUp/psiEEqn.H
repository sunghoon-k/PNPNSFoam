netCharge = zPlus*cPlus + zMinus*cMinus;
//Info << netCharge.dimensions() << endl;
fvScalarMatrix psiEEqn
(
    fvm::laplacian(psiE)
  //+ fvc::laplacian(psiE)
  //- fvc::surfaceIntegrate(gradpsiEf*magSf)
    //==
    //- netCharge/2/sqr(DebL_nonD)/sqr(l0_one)//*oneOverSquareDebL_nonD/sqr(l0_one)
);
psiEEqn.relax();

fvScalarMatrix psiIEqn
(
    fvm::laplacian(psiI)
);

psiIEqn.relax();



/*
if (currentIter != nOuterIter -1)
{
    psiEEqn.relax();
}
*/