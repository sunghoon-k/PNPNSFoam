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

    Info << "Reading field psiI\n" << endl;
    volScalarField psiI
    (
        IOobject
        (
            "psiI",
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
        fvc::snGrad(psiE) * mesh.magSf()
    );

    Info << "Reading field DPlusphiE\n" << endl;
    surfaceScalarField DPlusphiE
    (
        IOobject
        (
            "DPlusphiE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        DPlus_nonD * zPlus * phiE
    );    

    Info << "Reading field DMinusphiE\n" << endl;
    surfaceScalarField DMinusphiE
    (
        IOobject
        (
            "DMinusphiE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        DMinus_nonD * zMinus * phiE
    );    

    Info << "Reading field phiDcPlus\n" << endl;
    surfaceScalarField phiDcPlus
    (
        IOobject
        (
            "phiDcPlus",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        DPlus_nonD*zPlus*fvc::interpolate(cPlus)
    );

    Info << "Reading field phiDcMinus\n" << endl;
    surfaceScalarField phiDcMinus
    (
        IOobject
        (
            "phiDcMinus",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        DMinus_nonD*zMinus*fvc::interpolate(cMinus)
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
*/
    /*
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
            IOobject::AUTO_WRITE
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
            IOobject::NO_WRITE
        ),
        (zPlus*cPlus + zMinus*cMinus)
    );

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
