#include "readBlockSolverControls.H"
fvBlockMatrix<vector4> UpEqn(Up);

volScalarField cV
(
    IOobject
    (
        "cV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 3, 0, 0, 0, 0, 0), 0.0)
);
scalarField& cV_ = cV.internalField();
forAll(cV_, k) {cV_[k] = mesh.V()[k];}

dimensionedScalar dimL("dimL",dimensionSet(0,1,0,0,0,0,0), 1);
dimensionedScalar dimU("dimU",dimensionSet(0,1,-1,0,0,0,0), 1);

for(int NewtonIter = 0; NewtonIter < nNewtonIterNS; NewtonIter++)
{
    volVectorField dU(U);
    volScalarField dp(p);

    p.storePrevIter();

    #include "UNewtonEqn.H"
    #include "pNewtonEqn.H"
    #include "couplingTermsNSNewton.H"
    #include "calcFNS.H"

    // Solve the block matrix
    maxResidual_NS = cmptMax(UpEqn.solve().initialResidual());

    // Retrieve solution
    UpEqn.retrieveSolution(0, dU.internalField());
    UpEqn.retrieveSolution(3, dp.internalField());

    scalar lengthUx(0);
    scalar lengthUy(0);
    scalar lengthUz(0);
    scalar lengthp(0);

    // residual-length
    forAll(dU.internalField(),i)
    { 
        lengthUx   += sqr(dU.internalField()[i](0));
        lengthUy   += sqr(dU.internalField()[i](1)); 
        lengthUz   += sqr(dU.internalField()[i](2));
    }
    forAll(dp.internalField(),i){ lengthp += sqr(dp.internalField()[i]); }

    Info <<" Iter # "<< NewtonIter  << " - Residual: " <<"[" 
    << "Ux:" << Foam::sqrt(lengthUx) << "  " 
    << "Uy:" << Foam::sqrt(lengthUy) << "  " 
    << "Uz:" << Foam::sqrt(lengthUz) << "  " 
    << "p:" << Foam::sqrt(lengthp) << "]"  << endl;

    U -= dU;
    p -= dp;

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;
    #include "continuityErrs.H"
    Info << "maxResidual in NS: " << maxResidual_NS << endl;
    if (maxResidual_NS < convergenceCriterionNS)
    {
        Info << "reached criterion!" << endl; 
        break; 
    }

    p.relax();
}
