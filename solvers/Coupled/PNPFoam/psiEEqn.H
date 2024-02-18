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
if (i != nOuterIter -1)
{
  psiEEqn.relax();
}

psiE *= psiEScaling;
PNPEqn.insertEquation(0, psiEEqn);
