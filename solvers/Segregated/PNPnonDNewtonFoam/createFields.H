//    #include "NormalizedParameters.H"

    Info << "Reading field psiE\n" << endl;
    volScalarField psiE
    (
        IOobject
        (
            "psiE",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Reading field cPlus\n" << endl;
    volScalarField cPlus
    (
        IOobject
        (
            "cPlus",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Reading field cMinus\n" << endl;
    volScalarField cMinus
    (
        IOobject
        (
            "cMinus",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
/*
    Info << "Reading field netCharge\n" << endl;
    volScalarField netCharge
    (
        IOobject
        (
            "netCharge",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        F*(zPlus*cPlus + zMinus*zMinus)
    );

*/


/**********************************************************************/
    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    ElectrochemicalSystem ECsystem(mesh, transportProperties);
    #include "NormalizedParameters.H"

    word phiDyMBound(transportProperties.subDict("phiDyMBoundary").lookup("name"));

    // dimensionedScalar mobilityPlus((e*DPlus*zPlus)/(kB*T));
    // dimensionedScalar mobilityMinus((e*DMinus*zMinus)/(kB*T));

    Info << "Reading field psiEGrad\n" << endl;
    volVectorField psiEGrad
    (
        IOobject
        (
            "psiEGrad",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(psiE)
    );

    Info << "Reading field phiE\n" << endl;
    surfaceScalarField phiE
    (
        IOobject
        (
            "phiE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::snGrad(psiE) * mesh.magSf() // dimensions() = [0 -1 0 0 0 0 0] + [0 2 0 0 0 0 0] = [0 1 0 0 0 0 0]
    );

    
    Info << "Reading field currentcPlus\n" << endl;
    volVectorField currentcPlus
    (
        IOobject
        (
            "currentcPlus",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -( DPlus_nonD * fvc::grad(cPlus) + DPlus_nonD * psiEGrad * cPlus)
    );
    
    Info << "Reading field currentcMinus\n" << endl;
    volVectorField currentcMinus
    (
        IOobject
        (
            "currentcMinus",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -( DMinus_nonD * fvc::grad(cMinus) + DMinus_nonD * psiEGrad * cMinus)
    );

    Info << "Reading field netCharge\n" << endl;
    volScalarField netCharge
    (
        IOobject
        (
            "netCharge",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (zPlus*cPlus + zMinus*cMinus)
    );
