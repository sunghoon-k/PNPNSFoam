// Pressure parts of the continuity equation
surfaceScalarField presSource
(
    "presSource",
    fvc::interpolate(rAU)*
    (
      (
        fvc::interpolate(fvc::grad(p)) 
      //+ 
      //  fvc::interpolate(netCharge) * fvc::interpolate(fvc::grad(psiE))
      )& mesh.Sf()
    )
);

fvScalarMatrix pEqn
(
  - fvm::laplacian(rAU, p)
 ==
  - fvc::div(presSource)
);

pEqn.setReference(pRefCell, pRefValue);
