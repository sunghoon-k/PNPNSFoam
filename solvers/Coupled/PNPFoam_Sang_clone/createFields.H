struct physicalVolumeProperties{
word name;
scalar fixedCharge;
scalar diffusivity;
scalar transportNumber;
};
struct ionProperties{
word name;
scalar chargeNumber;
scalar diffusivity;

};
struct Ion{
word name;
scalar z;
scalar D;
};


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

Info<< "Reading field Phi\n" << endl;
volScalarField Phi
(
    IOobject
    (
        "Phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
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

Info<< "Reading field p\n" << endl;
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

Info << "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(U) & mesh.Sf()
);

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);


Info << "Creating field Up\n" << endl;
volVector4Field Up
(
    IOobject
    (
        "Up",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector4("zero", dimless, vector4::zero)
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
///***
Info<< "Creating field dC1\n" << endl;
volScalarField dC1
(
    IOobject
    (
        "dC1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("F1", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0)
);

Info<< "Creating field dC2\n" << endl;
volScalarField dC2
(
    IOobject
    (
        "dC2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("F2", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0)
);

Info<< "Creating field dPhi\n" << endl;
volScalarField dPhi
(
    IOobject
    (
        "dPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("F3", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0)
);

///***
Info<< "Reading field C1\n" << endl;
volScalarField C1
(
    IOobject
    (
        "C1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field C2\n" << endl;
volScalarField C2
(
    IOobject
    (
        "C2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);




Info<< "Creating field blockPNP\n" << endl;
volVector3Field blockPNP
(
    IOobject
    (
        "blockPNP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector3("zero",dimless, vector3::zero)
);

word phiDyMBound(transportProperties.subDict("phiDyMBoundary").lookup("name"));
