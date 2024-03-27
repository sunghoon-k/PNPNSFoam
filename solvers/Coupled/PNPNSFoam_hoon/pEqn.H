
fvScalarMatrix pEqn
(
  - fvm::laplacian(rAU, p)
// ==
//  - fvc::div(presSource)
);

pEqn.setReference(pRefCell, pRefValue);
