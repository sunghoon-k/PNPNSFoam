fvScalarMatrix psiEEqn
(
    fvm::laplacian(psiE)
  // + netCharge/(eps0*epsr)
  ==
  - fvc::laplacian(psiE_old)
);
if (currentIter != nOuterIter -1)
{
    psiEEqn.relax();
}
PNPEqn.insertEquation(0, psiEEqn);

