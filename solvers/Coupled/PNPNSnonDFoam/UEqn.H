volScalarField divPhi
(
    "divPhi",
    fvc::div(phi)
);

dimensionedScalar nu_one("nu_one", nu.dimensions(), 1.0);
dimensionedScalar bodyForce_one("bodyForce_one", dimensionSet(0,2,-2,0,0,0,0), 1.0);

netCharge = zPlus*cPlus + zMinus*cMinus;
if(splitpsi)
{
  bodyForce = - netCharge * fvc::grad(psiI) * oneOver2SquareDebL_nonD;
}
else
{
  bodyForce = - netCharge * fvc::grad(psiE) * oneOver2SquareDebL_nonD;
}
// Momentum equation
fvVectorMatrix UEqn
(
//        fvm::ddt(U)
    Re * fvm::div(phi, U)
  ==
    fvm::laplacian(nu_one, U)
  + bodyForce*bodyForce_one
);

if(solveTransient)
{
    UEqn += (1/Sc)*fvm::ddt(U);
}

rAU = 1.0/UEqn.A();

// Under-relax momentum. Note this will destroy the H and A
UEqn.relax();

