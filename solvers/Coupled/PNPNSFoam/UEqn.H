{    
    volScalarField divPhi
    (
        "divPhi",
        fvc::div(phi)
    );

    //volScalarField rhoEOverrho(netCharge/rho);
    // Info << "rhoEOverrho" << rhoEOverrho.dimensions() << endl;
    // Momentum equation
    fvVectorMatrix UEqn
    (
//        fvm::ddt(U)
      fvm::div(phi, U)
      ==
      fvm::laplacian(nu, U)
    //- rhoEOverrho * fvc::grad(psiE)
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
