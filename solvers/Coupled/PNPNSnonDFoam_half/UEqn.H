{
    // phi /= U0;
    
    volScalarField divPhi
    (
        "divPhi",
        fvc::div(phi)
    );

    dimensionedScalar ElectricBodyForce_one("ElectricBodyForce_one", (F*c0*psiE0/rho).dimensions(), 1.0);
    dimensionedScalar nu_one("nu_one", nu.dimensions(), 1.0);
    // Momentum equation
    fvVectorMatrix UEqn
    (
//        fvm::ddt(U)
        Re * fvm::div(phi, U)
      ==
        fvm::laplacian(nu_one, U)
      - ElectricBodyForce_one * (netCharge/sqr(DebL_nonD)) * fvc::grad(psiE)
    );

    if(solveTransient)
    {
        UEqn += (1/Sc)*fvm::ddt(U);
    }

    rAU = 1.0/UEqn.A();

    // Under-relax momentum. Note this will destroy the H and A
    UEqn.relax(1.0);

    NSEqn.insertEquation(0, UEqn);
}
