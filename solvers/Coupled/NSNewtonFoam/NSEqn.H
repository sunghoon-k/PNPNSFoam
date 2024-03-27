Info<< "########################################### "
    << "NSEqn"
    << " ###########################################" << nl << endl;            

p.storePrevIter();

volVectorField UT = U;
volScalarField pT = p;
surfaceScalarField presSource
(
    "presSource",
    fvc::interpolate(rAU)*
    (
      (
        fvc::interpolate(fvc::grad(p)) 
      //+ 
      //  fvc::interpolate(netCharge) * fvc::interpolate(fvc::grad(psiE))
      )& mesh.Sf()
    )
);

for(int NewtonIter = 0; NewtonIter < nNewtonIterNS; NewtonIter++)
{
    dimensionedScalar nu_one("nu_one", nu.dimensions(), 1.0);

    //volScalarField rAU_old = rAU;

    fvBlockMatrix<vector4> UpEqn(Up);

    #include "UEqn.H"
    #include "pEqn.H"
    #include "couplingTerms.H"
    presSource = fvc::interpolate(rAU)*((fvc::interpolate(fvc::grad(p)))& mesh.Sf());

    volVectorField FU = Re * fvc::div(phi, U) - fvc::laplacian(nu_one, U) + fvc::grad(p);
    volScalarField Fp = fvc::div(fvc::interpolate(U)&mesh.Sf()) - (fvc::laplacian(rAU, p) - fvc::div(presSource));

    if(pseudoTransient)
    {
        volVectorField FUt = (1/Sc) * fvc::ddt(U);
        FU += FUt;
    }

    // Assign F to blockM's right-hand side
    forAll(UpEqn.source(),cellI)
    {
        scalar cellVol = mesh.V()[cellI];
        UpEqn.source()[cellI](0) = -FU[cellI](0)*cellVol;
        UpEqn.source()[cellI](1) = -FU[cellI](1)*cellVol;
        UpEqn.source()[cellI](2) = -FU[cellI](2)*cellVol;
        UpEqn.source()[cellI](3) = -Fp[cellI]*cellVol;
    }

    maxResidual_NS = cmptMax(UpEqn.solve().initialResidual()); // residual = b - A*x_n

    // Retrieve solution
    UpEqn.retrieveSolution(0, U.internalField());
    UpEqn.retrieveSolution(3, p.internalField());

    scalar lengthUx(0);
    scalar lengthUy(0);
    scalar lengthUz(0);

    // residual-length
    forAll(U.internalField(),cellI)
    { 
        lengthUx += sqr(U.internalField()[cellI](0));
        lengthUy += sqr(U.internalField()[cellI](1)); 
        lengthUz += sqr(U.internalField()[cellI](2)); 
    }

    Info <<" Iter # "<< NewtonIter+1 << " - Residual: " <<"[" << "Ux:" << Foam::sqrt(lengthUx)
    << "  " << "Uy:" << Foam::sqrt(lengthUy) << "  " << "Uz:" << Foam::sqrt(lengthUz) << "]"  << endl;

    UT += U;
    pT += p;

    U = UT;
    p = pT;

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    //surfaceScalarField presSourcePhi = fvc::interpolate(rAU)*((fvc::interpolate(fvc::grad(p)))& mesh.Sf());
    presSource = fvc::interpolate(rAU)*((fvc::interpolate(fvc::grad(p)))& mesh.Sf());

    phi = (fvc::interpolate(U) & mesh.Sf()) - fvc::interpolate(rAU)*fvc::snGrad(p)*mesh.magSf() + presSource; // pEqn.flux() + presSource;

    #include "continuityErrs.H"
    
    scalar residual;
    residual = Foam::sqrt(lengthUx + lengthUy + lengthUz);
    Info << "residual: " << residual << endl;
    if (residual < convergenceCriterionNS)
    {
        Info << "break!!!!!!" << endl; 
        break; 
    }

    p.relax();
}
