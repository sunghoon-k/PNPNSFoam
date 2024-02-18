phiE = fvc::snGrad(psiE) * mesh.magSf();

fvScalarMatrix cPlusEqn
(
//    fvm::ddt(cPlus)
  - fvm::laplacian(DPlus, cPlus)
  - fvm::div(mobilityPlus * phiE, cPlus)
);

if(solveTransient)
{
  cPlusEqn += fvm::ddt(cPlus);
}

fvScalarMatrix cMinusEqn
(
//    fvm::ddt(cMinus)
  - fvm::laplacian(DMinus, cMinus)
  - fvm::div(mobilityMinus * phiE, cMinus)
);

if(solveTransient)
{
  cMinusEqn += fvm::ddt(cMinus);
}

if (i == nOuterIter -1)
{
  cPlusEqn.relax(1);
  cMinusEqn.relax(1);
}
else
{
  cPlusEqn.relax();
  cMinusEqn.relax();
}

PNPEqn.insertEquation(1, cPlusEqn);
PNPEqn.insertEquation(2, cMinusEqn);

//cPlusResidual = cPlusEqn.solve().initialResidual();
//cPlus.correctBoundaryConditions();

//cMinusResidual = cMinusEqn.solve().initialResidual();
//cMinus.correctBoundaryConditions();

//cMaxResidual = max(cPlusResidual, cMinusResidual);
