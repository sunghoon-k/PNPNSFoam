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

    word phiDyMBound(transportProperties.subDict("phiDyMBoundary").lookup("name"));
