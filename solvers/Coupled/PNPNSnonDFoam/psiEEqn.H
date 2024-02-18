//netCharge = (1/NA)*F*(zPlus*cPlus + zMinus*cMinus);
//Info << netCharge.dimensions() << endl;
fvScalarMatrix psiEEqn
(
    fvm::laplacian(psiE)
    //==
    //- netCharge/(eps0*epsr)
);
PNPNSEqn.insertEquation(4, psiEEqn);

/*
if (currentIter != nOuterIter -1)
{
    psiEEqn.relax();
}
*/