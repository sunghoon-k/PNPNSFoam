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

    Info << "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

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

    dimensionedScalar mobilityPlus((e*DPlus*zPlus)/(kB*T));
    dimensionedScalar mobilityMinus((e*DMinus*zMinus)/(kB*T));

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

/*
    Info << "Reading field phiMobility\n" << endl;
    surfaceScalarField phiMobilityPlus
    (
        IOobject
        (
            "phiMobilityPlus",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiE*(DPlus*e*zPlus)/(kB*T)
    );

    surfaceScalarField phiMobilityMinus
    (
        IOobject
        (
            "phiMobilityMinus",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiE*(DMinus*e*zMinus)/(kB*T)
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
            IOobject::NO_WRITE
        ),
        -( fvc::grad(cPlus) * DPlus + mobilityPlus * psiEGrad * cPlus)
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
            IOobject::NO_WRITE
        ),
        -( fvc::grad(cMinus) * DMinus + mobilityMinus* psiEGrad * cMinus)
    );

*/
    
    Info << "Reading field netCharge\n" << endl;
    volScalarField netCharge
    (
        IOobject
        (
            "netCharge",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (zPlus*cPlus + zMinus*cMinus)
    );

    Info << "Reading field bodyForce\n" << endl;
    volVectorField bodyForce
    (
        IOobject
        (
            "bodyForce",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        - netCharge * fvc::grad(psiE) / (sqr(DebL_nonD.value())*2)
    );

    Info<< "Creating field rAU\n" << endl;
    volScalarField rAU
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        runTime.deltaT()
    );

    mesh.schemesDict().setFluxRequired(p.name());
    mesh.schemesDict().setFluxRequired("pcorr");

    Info << "Creating field PNP\n" << endl;
    volVector3Field PNP
    (
        IOobject
        (
            "PNP",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector3("zero", dimless, vector3::zero)
    );

    Info << "Creating field NS\n" << endl;
    volVector4Field NS
    (
        IOobject
        (
            "PNPNS",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector4("zero", dimless, vector4::zero)
    );
