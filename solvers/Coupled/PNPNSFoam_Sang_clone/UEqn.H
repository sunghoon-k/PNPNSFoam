    volScalarField divPhi
    (
        "divPhi",
        fvc::div(phi)
    );

    scalar iEpsilon2 = 1.0/Epsilon2;
    dimensionedScalar one("one", dimensionSet(0, 2, -2, 0, 0, 0, 0), 1.0);

    fvVectorMatrix UEqn
    (
        Re*fvm::div(phi, U)
      - fvm::laplacian(nu, U)
      ==
        bodyForce*one//+ iEpsilon2*(z1*C1 + z2*C2)*fvc::grad(Phi)*one/2
    );

    rAU = 1.0/UEqn.A();

    UEqn.relax();

    UpEqn.insertEquation(0, UEqn);
