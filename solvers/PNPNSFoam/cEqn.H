/*
volScalarField rhoE(phiE * dimensionedScalar("norm1",dimless/dimArea,0.));       

rhoE += ( zPlus * cPlus * const_F );
rhoE += ( zMinus * cMinus * const_F );
*/

surfaceScalarField phiE(linearInterpolate(fvc::grad(psiE)) & mesh.Sf());

surfaceScalarField mobilityPlus((DPlus * e * zPlus/(kB * T))*phiE);
surfaceScalarField mobilityMinus((DMinus * e * zMinus/(kB * T))*phiE);

fvScalarMatrix cPlusEqn
(
//    fvm::ddt(cPlus)
    fvm::div(phi, cPlus)    
    ==
    fvm::laplacian(DPlus, cPlus)
  + fvm::div(mobilityPlus, cPlus)
);

if(solveTransient)
{
  cPlusEqn += fvm::ddt(cPlus);
}

fvScalarMatrix cMinusEqn
(
//    fvm::ddt(cMinus)
    fvm::div(phi, cMinus)    
    ==
    fvm::laplacian(DMinus, cMinus)
  + fvm::div(mobilityMinus, cMinus)
);

if(solveTransient)
{
  cMinusEqn += fvm::ddt(cMinus);
}

cPlusEqn.relax();
cMinusEqn.relax();

cPlusEqn *= cScaling;
cMinusEqn *= cScaling;

PNPNSEqn.insertEquation(5, cPlusEqn);
PNPNSEqn.insertEquation(6, cMinusEqn);