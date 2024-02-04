/*
volScalarField rhoE(phiE * dimensionedScalar("norm1",dimless/dimArea,0.));       

rhoE += ( zPlus * cPlus * const_F );
rhoE += ( zMinus * cMinus * const_F );
*/

fvScalarMatrix phiEEqn
(
    fvm::laplacian(phiE)
/* ==
  - rhoE/(epsilon0*epsilonr)*/
);

phiE_c_UpEqn.insertEquation(4, phiEEqn);
