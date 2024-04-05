volScalarField divPhi
(
    "divPhi",
    fvc::div(phi)
);


netCharge = zPlus*cPlus + zMinus*cMinus;

bodyForce = - netCharge * fvc::grad(psiE) * oneOver2SquareDebL_nonD;

// Momentum equation
fvVectorMatrix UEqn
(
//        fvm::ddt(U)
    Re * fvm::div(phi, U)
  ==
    fvm::laplacian(nu_one, U)
  + bodyForce*bodyForce_one
);

if(pseudoTransient)
{
    UEqn += (1/Sc)*fvm::ddt(U);
}

rAU = 1.0/UEqn.A();

// Under-relax momentum. Note this will destroy the H and A
UEqn.relax();

