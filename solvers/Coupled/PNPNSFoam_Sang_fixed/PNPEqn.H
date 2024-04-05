volScalarField cPlusT = cPlus;
volScalarField cMinusT = cMinus;
volScalarField psiET = psiE;

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

for(int NewtonIter = 0; NewtonIter < nNewtonIterPNP; NewtonIter++)
{
//    volScalarField dcPlus(cPlus);
//    volScalarField dcMinus(cMinus);
//    volScalarField dpsiE(psiE);    
    
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
    #include "couplingTermsPNP.H"

    // Calculate Right-Hand Side F
    volScalarField zPluscPlusVol = zPlus*cPlus*oneOver2SquareDebL_nonD; // *cV
    volScalarField zMinuscMinusVol = zMinus*cMinus*oneOver2SquareDebL_nonD; // *cV
/*
    volScalarField F1 = (fvc::laplacian(-DPlus_nonD, cPlus) + fvc::div(DPluszPlusphiE, cPlus))*cV; //  + fvc::laplacian(DPluszPluscPlusf, psiE)
    volScalarField F2 = (fvc::laplacian(-DMinus_nonD, cMinus) + fvc::div(DMinuszMinusphiE, cMinus) )*cV; // + fvc::laplacian(DMinuszMinuscMinusf, psiE)
    volScalarField F3 = (fvc::laplacian(psiE) + (zPlus*cPlus + zMinus*cMinus)*oneOver2SquareDebL_nonD/sqr(dimL))*cV;

*/
    volScalarField F1 = fvc::surfaceIntegrate(jPlus)*cV;// *sqr(dimL);
    volScalarField F2 = fvc::surfaceIntegrate(jMinus)*cV;// *sqr(dimL);
    volScalarField F3 = fvc::surfaceIntegrate(gradpsiEf*magSf)*cV + (zPluscPlusVol + zMinuscMinusVol)/sqr(dimL)*cV; // *cV


    // Assign F to blockM's right-hand side
    forAll(PNPEqn.source(),cellI)
    {
        PNPEqn.source()[cellI](0) = -F1[cellI];
        PNPEqn.source()[cellI](1) = -F2[cellI];
        PNPEqn.source()[cellI](2) = -F3[cellI];
    }

    maxResidual_PNP = cmptMax(PNPEqn.solve().initialResidual()); // residual = b - A*x_n

    // Retrieve solution
    PNPEqn.retrieveSolution(0, cPlus.internalField());
    PNPEqn.retrieveSolution(1, cMinus.internalField());
    PNPEqn.retrieveSolution(2, psiE.internalField());
    //-Calculate and screen residual-length
    // Initialize residual-length

    scalar lengthcPlus(0);
    scalar lengthcMinus(0);
    scalar lengthpsiE(0);

    // residual-length
    forAll(cPlus.internalField(),i){ lengthcPlus   += sqr(cPlus.internalField()[i]); }
    forAll(cMinus.internalField(),i){ lengthcMinus += sqr(cMinus.internalField()[i]); }
    forAll(psiE.internalField(),i){ lengthpsiE     += sqr(psiE.internalField()[i]); }

    Info <<" Iter # "<< NewtonIter  << " - Residual: " <<"[" << "cPlus:" << Foam::sqrt(lengthcPlus)
    << "  " << "cMinus:" << Foam::sqrt(lengthcMinus) << "  " << "psiE:" << Foam::sqrt(lengthpsiE) << "]"  << endl;

// Correct C1,C2,Phi
    cPlusT   += cPlus;
    cMinusT  += cMinus;
    psiET    += psiE;

    // Assign C1T,C2T,PhiT to C1,C2,Phi
    cPlus  = cPlusT;
    cMinus  = cMinusT;
    psiE = psiET;
    
    /// CorrectBoundaryCondition
    cPlus.correctBoundaryConditions();
    cMinus.correctBoundaryConditions();
    psiE.correctBoundaryConditions();
    // Correct Boundary value according to J-Flux

    //#include "boundC.H"
    scalar residual;
    residual = Foam::sqrt(lengthcPlus + lengthcMinus + lengthpsiE);
    Info << "residual: " << residual << endl;
    if (residual < convergenceCriterionPNP)
    {
        Info << "break!!!!!!" << endl; 
        break; 
    }
}

