{
    // Bound the pressure
    dimensionedScalar p1 = min(p);
    dimensionedScalar p2 = max(p);

    if (p1 < pMin)
    {
        Info<< "p: " << p1.value() 
            << ".  Bounding." << endl;

        p.max(pMin); // max(p, pMin)
        p.correctBoundaryConditions();
    }

/*
    // Bound the velocity
    volScalarField magU = mag(U);
    dimensionedScalar U1 = max(magU);

    if (U1 > UMax)
    {
        Info<< "U: " << U1.value() << ".  Bounding." << endl;

        // pos: return 1, if input >= 0
        // pos(magU - UMax) == 1, if |U| >= UMax
        // neg: return 1, if input < 0
        // neg(magU - UMax) == 0, if |U| >= UMax
        // Ulimiter = UMax/(|U|+smallU) or 1
        volScalarField Ulimiter = pos(magU - UMax)*UMax/(magU + smallU)
            + neg(magU - UMax);
        Ulimiter.max(scalar(0));
        Ulimiter.min(scalar(1));
        // 0 <= Ulimiter <= 1

        U *= Ulimiter;
        U.correctBoundaryConditions();
    }

*/

}
