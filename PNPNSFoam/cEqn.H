/*
volScalarField rhoE(phiE * dimensionedScalar("norm1",dimless/dimArea,0.));       

rhoE += ( zPlus * cPlus * const_F );
rhoE += ( zMinus * cMinus * const_F );
*/

surfaceScalarField phiE(linearInterpolate(fvc::grad(psiE)) & mesh.Sf());

surfaceScalarField mobilityPlus((DPlus*e_const*zPlus/(k_const*T_const))*phiE);
surfaceScalarField mobilityMinus((DMinus*e_const*zMinus/(k_const*T_const))*phiE);

fvScalarMatrix cPlusEqn
(
    fvm::ddt(cPlus)
  + fvm::div(phi, cPlus)    
    ==
    fvm::laplacian(DPlus, cPlus)
  + fvm::div(mobilityPlus, cPlus)
);

fvScalarMatrix cMinusEqn
(
    fvm::ddt(cMinus)
  + fvm::div(phi, cMinus)    
    ==
    fvm::laplacian(DMinus, cMinus)
  + fvm::div(mobilityMinus, cMinus)
);

cPlusEqn.relax();
cMinusEqn.relax();

PNPNSEqn.insertEquation(5, cPlusEqn);
PNPNSEqn.insertEquation(6, cMinusEqn);