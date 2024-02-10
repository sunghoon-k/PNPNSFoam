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

#   include "createPhi.H"
// #   include "createPhiE.H"

/*
    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );
*/
    // Block vector field for velocity (first entry) and pressure (second
    // entry).

    Info << "Creating field PNPNS\n" << endl;
    volVector7Field PNPNS
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
        dimensionedVector7("zero", dimless, vector7::zero)
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

/*
    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    dimensionedScalar rho
    (
        transportProperties.lookup("rho")
    );

    dimensionedScalar zPlus
    (
        transportProperties.lookup("zPlus")
    );

    dimensionedScalar zMinus
    (
        transportProperties.lookup("zMinus")
    );

    dimensionedScalar DPlus
    (
        transportProperties.lookup("DPlus")
    );

    dimensionedScalar DMinus
    (
        transportProperties.lookup("DMinus")
    );

    dimensionedScalar epsilon0
    (
        transportProperties.lookup("epsilon0")
    );

    dimensionedScalar epsilonr
    (
        transportProperties.lookup("epsilonr")
    );

    dimensionedScalar F_const // Faraday's constant
    (
        transportProperties.lookup("F_const")
    );

    dimensionedScalar e_const // elementary charge
    (
        transportProperties.lookup("e_const")
    );

    dimensionedScalar k_const // Boltzmann's constant
    (
        transportProperties.lookup("k_const")
    );

    dimensionedScalar T_const // Boltzmann's constant
    (
        transportProperties.lookup("T_const")
    );

*/

    word phiDyMBound(transportProperties.subDict("phiDyMBoundary").lookup("name"));
