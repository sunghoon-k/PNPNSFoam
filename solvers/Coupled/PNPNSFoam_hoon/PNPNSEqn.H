volScalarField cPlusT = cPlus;
volScalarField cMinusT = cMinus;
volScalarField psiET = psiE;
volVectorField UT = U;
volScalarField pT = p;

dimensionedScalar dimL("dimL",dimensionSet(0,1,0,0,0,0,0), 1);
dimensionedScalar dimU("dimU",dimensionSet(0,1,-1,0,0,0,0), 1);

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

for(int NewtonIter = 0; NewtonIter < nNewtonIter; NewtonIter++)
{
    p.storePrevIter();


    /********************************************************************************/

    surfaceScalarField gradcPlusf   = fvc::snGrad(cPlus);            
    surfaceScalarField gradcMinusf  = fvc::snGrad(cMinus);            
    surfaceScalarField cPlusf       = fvc::interpolate(cPlus);            
    surfaceScalarField cMinusf      = fvc::interpolate(cMinus);            
    surfaceScalarField gradpsiEf    = fvc::snGrad(psiE);    
    surfaceScalarField gradpsiIf    = fvc::snGrad(psiI);    
    surfaceScalarField magSf        = mag(mesh.Sf());            

    // Calculate J-Flux
    surfaceScalarField jPlus = -DPlus_nonD*(gradcPlusf + zPlus*cPlusf*gradpsiEf)*magSf      + Pe*cPlusf*phi/(dimL*dimU);
    surfaceScalarField jMinus = -DMinus_nonD*(gradcMinusf + zMinus*cMinusf*gradpsiEf)*magSf + Pe*cMinusf*phi/(dimL*dimU);

    // Assign J-Flux to specific value
    forAll(cPlus.boundaryField(),patchI)
    {
        if(cPlus.boundaryField()[patchI].type() == "fixedIonicFlux" or cPlus.boundaryField()[patchI].type() == "zeroIonicFlux_nonD")
        {
            scalarField& jPlusf = jPlus.boundaryField()[patchI];
            forAll(jPlusf, i){ jPlusf[i] = 0; }
        }
    }

    forAll(cMinus.boundaryField(),patchI)
    {
        if(cMinus.boundaryField()[patchI].type() == "fixedIonicFlux" or cMinus.boundaryField()[patchI].type() == "zeroIonicFlux_nonD")
        {
            scalarField& jMinusf = jMinus.boundaryField()[patchI];
            forAll(jMinusf, i){ jMinusf[i] = 0; }
        }
    }

    /// Jacobian-Matrix Assembling
    surfaceScalarField DPluszPlusphiE("DPluszPlusphiE",             -DPlus_nonD*zPlus*(gradpsiEf+gradpsiIf)*magSf   + Pe*cPlusf*phi/(dimU*l0_one)); // 
    surfaceScalarField DMinuszMinusphiE("DMinuszMinusphiE",         -DMinus_nonD*zMinus*(gradpsiEf+gradpsiIf)*magSf + Pe*cMinusf*phi/(dimU*l0_one));
    surfaceScalarField DPluszPluscPlusf("DPluszPluscPlusf",         -DPlus_nonD*zPlus*cPlusf);
    surfaceScalarField DMinuszMinuscMinusf("DMinuszMinuscMinusf",   -DMinus_nonD*zMinus*cMinusf);


    forAll(cPlus.boundaryField(),patchI)
    {
        if(cPlus.boundaryField()[patchI].type() == "fixedIonicFlux" or cPlus.boundaryField()[patchI].type() == "zeroIonicFlux_nonD")
        {
            scalarField& DPluszPlusphiE_ = DPluszPlusphiE.boundaryField()[patchI];
            forAll(DPluszPlusphiE_, i) { DPluszPlusphiE_[i] = 0; }
            scalarField& DPluszPluscPlusf_ = DPluszPluscPlusf.boundaryField()[patchI];
            forAll(DPluszPluscPlusf_, i) { DPluszPluscPlusf_[i] = 0; }
        }
    }

    forAll(cMinus.boundaryField(),patchI)
    {
        if(cMinus.boundaryField()[patchI].type() == "fixedIonicFlux" or cMinus.boundaryField()[patchI].type() == "zeroIonicFlux_nonD")
        {
            //Info << "zeroIonicFlux_nonD ############################"<<endl;
            scalarField& DMinuszMinusphiE_ = DMinuszMinusphiE.boundaryField()[patchI];
            forAll(DMinuszMinusphiE_, i) { DMinuszMinusphiE_[i] = 0; }
            scalarField& DMinuszMinuscMinusf_ = DMinuszMinuscMinusf.boundaryField()[patchI];
            forAll(DMinuszMinuscMinusf_, i) { DMinuszMinuscMinusf_[i] = 0; }
        }
    }

    #include "cEqn.H"
    #include "psiEEqn.H"
    #include "UEqn.H"
    #include "pEqn.H"
    #include "couplingTermsPNPNS.H"

    // Calculate Right-Hand Side F
    volScalarField zPluscPlusVol = zPlus*cPlus*oneOver2SquareDebL_nonD; // *cV
    volScalarField zMinuscMinusVol = zMinus*cMinus*oneOver2SquareDebL_nonD; // *cV

    volScalarField F1 = fvc::surfaceIntegrate(jPlus);// *sqr(dimL);
    volScalarField F2 = fvc::surfaceIntegrate(jMinus);// *sqr(dimL);
    volScalarField F3 = fvc::surfaceIntegrate(gradpsiEf*magSf) + (zPluscPlusVol + zMinuscMinusVol)/sqr(dimL); // *cV
    volVectorField FU = Re * fvc::div(phi, U) - fvc::laplacian(nu_one, U) + fvc::grad(p);
    
    presSource = fvc::interpolate(rAU)*((fvc::interpolate(fvc::grad(p)))& mesh.Sf());
    volScalarField Fp = fvc::div(fvc::interpolate(U)&mesh.Sf()) - (fvc::laplacian(rAU, p) - fvc::div(presSource));

    if(pseudoTransient)
    {
        volScalarField F1t = ddtc_one * fvc::ddt(cPlus);
        volScalarField F2t = ddtc_one * fvc::ddt(cMinus);
        F1 += F1t;
        F2 += F2t;

        volVectorField FUt = (1/Sc) * fvc::ddt(U);
        FU += FUt;
    }

    if(solveTransient)
    {
        volScalarField F1t = ddtc_one * ddtcPlus;
        volScalarField F2t = ddtc_one * ddtcMinus;
        F1 += F1t;
        F2 += F2t;

        volVectorField FUt = (1/Sc) * ddtU;
        FU += FUt;
    }

    // Assign F to blockM's right-hand side
    forAll(PNPNSEqn.source(),cellI)
    {
        scalar cellVol = mesh.V()[cellI];

        PNPNSEqn.source()[cellI](0) = -F1[cellI]*cellVol;
        PNPNSEqn.source()[cellI](1) = -F2[cellI]*cellVol;
        PNPNSEqn.source()[cellI](2) = -F3[cellI]*cellVol;
        PNPNSEqn.source()[cellI](3) = -FU[cellI](0)*cellVol;
        PNPNSEqn.source()[cellI](4) = -FU[cellI](1)*cellVol;
        PNPNSEqn.source()[cellI](5) = -FU[cellI](2)*cellVol;
        PNPNSEqn.source()[cellI](6) = -Fp[cellI]*cellVol;
    }

    maxResidual = cmptMax(PNPNSEqn.solve().initialResidual()); // residual = b - A*x_n

    // Retrieve solution
    PNPNSEqn.retrieveSolution(0, cPlus.internalField());
    PNPNSEqn.retrieveSolution(1, cMinus.internalField());
    PNPNSEqn.retrieveSolution(2, psiE.internalField());
    PNPNSEqn.retrieveSolution(3, U.internalField());
    PNPNSEqn.retrieveSolution(6, p.internalField());
    //-Calculate and screen residual-length
    // Initialize residual-length

    scalar lengthcPlus(0);
    scalar lengthcMinus(0);
    scalar lengthpsiE(0);
    scalar lengthUx(0);
    scalar lengthUy(0);
    scalar lengthUz(0);

    // residual-length
    forAll(cPlus.internalField(),i){ lengthcPlus   += sqr(cPlus.internalField()[i]); }
    forAll(cMinus.internalField(),i){ lengthcMinus += sqr(cMinus.internalField()[i]); }
    forAll(psiE.internalField(),i){ lengthpsiE     += sqr(psiE.internalField()[i]); }
    forAll(U.internalField(),cellI)
    { 
        lengthUx += sqr(U.internalField()[cellI](0));
        lengthUy += sqr(U.internalField()[cellI](1)); 
        lengthUz += sqr(U.internalField()[cellI](2)); 
    }

    Info <<" Iter # "<< NewtonIter  << " - Residual: " <<"[" 
    << "cPlus:" << Foam::sqrt(lengthcPlus) << "  " 
    << "cMinus:" << Foam::sqrt(lengthcMinus) << "  " 
    << "psiE:" << Foam::sqrt(lengthpsiE) << "  "
    << "Ux:" << Foam::sqrt(lengthUx) << "  " 
    << "Uy:" << Foam::sqrt(lengthUy) << "  " 
    << "Uz:" << Foam::sqrt(lengthUz) 
    << "]"  << endl;

// Correct C1,C2,Phi
    cPlusT  += cPlus;
    cMinusT += cMinus;
    psiET   += psiE;
    UT      += U;
    pT      += p;

    // Assign C1T,C2T,PhiT to C1,C2,Phi
    cPlus   = cPlusT;
    cMinus  = cMinusT;
    psiE    = psiET;
    U       = UT;
    p       = pT;

    /// CorrectBoundaryCondition
    cPlus.correctBoundaryConditions();
    cMinus.correctBoundaryConditions();
    psiE.correctBoundaryConditions();
    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    presSource = fvc::interpolate(rAU)*((fvc::interpolate(fvc::grad(p)))& mesh.Sf());

    phi = (fvc::interpolate(U) & mesh.Sf()) - fvc::interpolate(rAU)*fvc::snGrad(p)*mesh.magSf() + presSource; // pEqn.flux() + presSource;

//    phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;
    #include "continuityErrs.H"

    //#include "boundC.H"
    scalar residual;
    residual = Foam::sqrt(lengthcPlus + lengthcMinus + lengthpsiE + lengthUx + lengthUy + lengthUz);
    Info << "PNPNS residual: " << residual << endl;
    if (residual < convergenceCriterion)
    {
        Info << "PNPNS residual < " << convergenceCriterion << nl 
        << "break!!!!!!" << endl; 
        break; 
    }
    p.relax();

}

