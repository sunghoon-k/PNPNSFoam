    Info << "Reading field phiE\n" << endl;
    volScalarField phiE
    (
        IOobject
        (
            "phiE",
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

/*
    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );
*/
    // Block vector field for velocity (first entry) and pressure (second
    // entry).

/*
    typedef Field<vector7> vector7Field;
    typedef Field<sphericalTensor7> sphericalTensor7Field;
    typedef Field<diagTensor7> diagTensor7Field;
    typedef Field<tensor7> tensor7Field;
*/

    typedef GeometricField<vector7, fvPatchField, volMesh> volVector7Field;
    typedef GeometricField<sphericalTensor7, fvPatchField, volMesh> volSphericalTensor7Field;
    typedef GeometricField<diagTensor7, fvPatchField, volMesh> volDiagTensor7Field;
    typedef GeometricField<tensor7, fvPatchField, volMesh> volTensor7Field;

    typedef dimensioned<vector7> dimensionedVector7;
    typedef dimensioned<sphericalTensor7> dimensionedSphericalTensor7;
    typedef dimensioned<diagTensor7> dimensionedDiagTensor7;
    typedef dimensioned<tensor7> dimensionedTensor7;

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

    dimensionedScalar F_const
    (
        transportProperties.lookup("const_F")
    );
