/*
volScalarField rhoE(phiE * dimensionedScalar("norm1",dimless/dimArea,0.));       

rhoE += ( zPlus * cPlus * const_F );
rhoE += ( zMinus * cMinus * const_F );
*/

fvScalarMatrix psiEEqn
(
    fvm::laplacian(psiE)
  /*
  ==
  - fvm::SuSp((F_const*zPlus)/(epsilon0*epsilonr), cPlus)
  - fvm::SuSp((F_const*zMinus)/(epsilon0*epsilonr), cMinus)

  */
);

psiEEqn.relax();
psiE *= psiEScaling;
PNPNSEqn.insertEquation(4, psiEEqn);
