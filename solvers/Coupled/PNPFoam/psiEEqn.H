//netCharge = e*(zPlus*cPlus + zMinus*cMinus);
//Info << netCharge.dimensions() << endl;
fvScalarMatrix psiEEqn
(
    fvm::laplacian(psiE)
    //==
    //- netCharge/(eps0*epsr)
);
psiEEqn.relax();
PNPEqn.insertEquation(2, psiEEqn);

scalarField PhicouplingPlus(mesh.nCells(),((e*zPlus)/(eps0*epsr)).value());
scalarField PhicouplingMinus(mesh.nCells(),((e*zMinus)/(eps0*epsr)).value());
PNPEqn.insertEquationCoupling(2, 0, PhicouplingPlus); // 0: psiE <- 1: cPlus
PNPEqn.insertEquationCoupling(2, 1, PhicouplingMinus); // 0: psiE <- 2: cMinus

/*
if (currentIter != nOuterIter -1)
{
    psiEEqn.relax();
}
*/