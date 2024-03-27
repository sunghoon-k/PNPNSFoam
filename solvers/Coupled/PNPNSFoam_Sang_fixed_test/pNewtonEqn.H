// Pressure parts of the continuity equation
surfaceScalarField dpresSource
(
    "dpresSource",
    fvc::interpolate(rAU)*
    (
      (
        fvc::interpolate(fvc::grad(dp)) 
      //+ 
      //  fvc::interpolate(netCharge) * fvc::interpolate(fvc::grad(psiE))
      )& mesh.Sf()
    )
);

fvScalarMatrix pEqn
(
  - fvm::laplacian(rAU, dp)
 ==
  - fvc::div(dpresSource)
);

pEqn.setReference(pRefCell, pRefValue);