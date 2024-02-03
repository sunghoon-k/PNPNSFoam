{
    volScalarField divPhi
    (
        "divPhi",
        fvc::div(phi)
    );

    // Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + fvm::laplacian(nu, U)
      // + fvm::grad(rhoE/rho, phiE)
    );

    rAU = 1.0/UEqn.A();

    // Under-relax momentum. Note this will destroy the H and A
    UEqn.relax();

    PNPNSEqn.insertEquation(0, UEqn);
}
