{
    volScalarField divPhi
    (
        "divPhi",
        fvc::div(phi)
    );

    dimensionedScalar nu_one("nu_one", nu.dimensions(), 1.0);
    // Momentum equation
    fvVectorMatrix UEqn
    (
//        fvm::ddt(U)
        Re * fvm::div(phi, U)
      ==
        fvm::laplacian(nu_one, U)
    );

    if(solveTransient)
    {
        UEqn += (Re/Pe)*fvm::ddt(U);
    }

    rAU = 1.0/UEqn.A();

    // Under-relax momentum. Note this will destroy the H and A
    UEqn.relax(1.0);

    NSEqn.insertEquation(0, UEqn);
}
