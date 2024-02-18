Info<< "Reading/calculating face flux field phiE\n" << endl;

surfaceScalarField phiE
(
    IOobject
    (
        "phiE",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(fvc::grad(psiE)) & mesh.Sf()
);


// ************************************************************************* //