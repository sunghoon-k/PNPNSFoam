// Pressure parts of the continuity equation
{

  fvScalarMatrix pEqn
  (
    - fvm::laplacian(rAU, p)
  //==
  //  - fvc::div(presSource)
  );

  pEqn.setReference(pRefCell, pRefValue);
  UpEqn.insertEquation(3, pEqn);
}