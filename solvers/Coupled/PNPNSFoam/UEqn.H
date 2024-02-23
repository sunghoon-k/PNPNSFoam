{
    // phi /= U0;
    
    volScalarField divPhi
    (
        "divPhi",
        fvc::div(phi)
    );

    // Momentum equation
    fvVectorMatrix UEqn
    (
//        fvm::ddt(U)
      fvm::div(phi, U)
      ==
      fvm::laplacian(nu, U)
      // + fvm::grad(rhoE/rho, phiE)
    );

    if(solveTransient)
    {
        UEqn += fvm::ddt(U);
    }

    rAU = 1.0/UEqn.A();

    // Under-relax momentum. Note this will destroy the H and A
    UEqn.relax(1.0);

    PNPNSEqn.insertEquation(0, UEqn);
}
