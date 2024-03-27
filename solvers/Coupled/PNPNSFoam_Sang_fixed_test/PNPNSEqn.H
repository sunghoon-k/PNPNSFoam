volScalarField cV( IOobject ( "cV", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ), mesh,
    dimensionedScalar("zero", dimensionSet(0, 3, 0, 0, 0, 0, 0), 0.0)
);
scalarField& cV_ = cV.internalField();
forAll(cV_, k) {cV_[k] = mesh.V()[k];}

dimensionedScalar dimL("dimL",dimensionSet(0,1,0,0,0,0,0), 1);
dimensionedScalar dimU("dimU",dimensionSet(0,1,-1,0,0,0,0), 1);
for(int NewtonIter = 0; NewtonIter < nNewtonIterPNP; NewtonIter++)
{
    /********************************************************************************/

    volScalarField dcPlus(cPlus);
    volScalarField dcMinus(cMinus);
    volScalarField dpsiE(psiE);    
    
    /********************************************************************************/

    surfaceScalarField gradcPlusf   = fvc::snGrad(cPlus);            
    surfaceScalarField gradcMinusf  = fvc::snGrad(cMinus);            
    surfaceScalarField cPlusf       = fvc::interpolate(cPlus);            
    surfaceScalarField cMinusf      = fvc::interpolate(cMinus);            
    surfaceScalarField gradpsiEf    = fvc::snGrad(psiE);    
    surfaceScalarField gradpsiIf    = fvc::snGrad(psiI);    
    surfaceScalarField magSf        = mag(mesh.Sf());            

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
    volScalarField F1 = (fvc::laplacian(-DPlus_nonD, cPlus) + fvc::div(DPluszPlusphiE, cPlus) + fvc::laplacian(DPluszPluscPlusf, psiE))*cV;
    volScalarField F2 = (fvc::laplacian(-DMinus_nonD, cMinus) + fvc::div(DMinuszMinusphiE, cMinus) + fvc::laplacian(DMinuszMinuscMinusf, psiE))*cV;
    volScalarField F3 = (fvc::laplacian(psiE) + (zPlus*cPlus + zMinus*cMinus)*oneOver2SquareDebL_nonD/sqr(dimL))*cV;

    if(debugMode)
    {
        // Before Assigning F to blockM's right-hand side
        forAll(PNPEqn.source(),cellI)
        {
            Info<<"cPlus.source: "<<PNPEqn.source()[cellI](0)<<endl;
            Info<<"cMinus.source: "<<PNPEqn.source()[cellI](1)<<endl;
            Info<<"psiE.source: "<<PNPEqn.source()[cellI](2)<<endl;

        }
    }

    // Assign F to blockM's right-hand side
    forAll(PNPEqn.source(),cellI)
    {
        PNPEqn.source()[cellI](0) += F1[cellI];
        PNPEqn.source()[cellI](1) += F2[cellI];
        PNPEqn.source()[cellI](2) += F3[cellI];
    }

    if(debugMode)
    {
        // Before Assigning F to blockM's right-hand side
        forAll(PNPEqn.source(),cellI)
        {
            Info<<"cPlus.source: "<<PNPEqn.source()[cellI](0)<<endl;
            Info<<"cMinus.source: "<<PNPEqn.source()[cellI](1)<<endl;
            Info<<"psiE.source: "<<PNPEqn.source()[cellI](2)<<endl;

        }
    }
    
    maxResidual_PNP = cmptMax(PNPEqn.solve().initialResidual()); // residual = b - A*x_n

    // Retrieve solution
    PNPEqn.retrieveSolution(0, dcPlus.internalField());
    PNPEqn.retrieveSolution(1, dcMinus.internalField());
    PNPEqn.retrieveSolution(2, dpsiE.internalField());
    //-Calculate and screen residual-length
    // Initialize residual-length

    scalar lengthcPlus(0);
    scalar lengthcMinus(0);
    scalar lengthpsiE(0);

    // residual-length
    forAll(dcPlus.internalField(),i){ lengthcPlus   += sqr(dcPlus.internalField()[i]); }
    forAll(dcMinus.internalField(),i){ lengthcMinus += sqr(dcMinus.internalField()[i]); }
    forAll(dpsiE.internalField(),i){ lengthpsiE     += sqr(dpsiE.internalField()[i]); }

    Info <<" Iter # "<< NewtonIter  << " - Residual: " <<"[" << "cPlus:" << Foam::sqrt(lengthcPlus)
    << "  " << "cMinus:" << Foam::sqrt(lengthcMinus) << "  " << "psiE:" << Foam::sqrt(lengthpsiE) << "]"  << endl;
/*
    

*/
    cPlus   -= dcPlus;
    cMinus  -= dcMinus;
    psiE    -= dpsiE;
    
    
/*
// Correct C1,C2,Phi
    cPlusT   += cPlus;
    cMinusT  += cMinus;
    psiET    += psiE;

    // Assign C1T,C2T,PhiT to C1,C2,Phi
    cPlus  = cPlusT;
    cMinus  = cMinusT;
    psiE = psiET;
*/

    /// CorrectBoundaryCondition
    cPlus.correctBoundaryConditions();
    cMinus.correctBoundaryConditions();
    psiE.correctBoundaryConditions();
    // Correct Boundary value according to J-Flux
    forAll(cPlus.boundaryField(),patchI)
    {
        if(cPlus.boundaryField()[patchI].type() == "fixedIonicFlux")
        {
            scalarField& Cb1 = cPlus.boundaryField()[patchI];
            const tmp<scalarField>&  Ci1 = cPlus.boundaryField()[patchI].patchInternalField();
            scalarField& Phib1 = psiE.boundaryField()[patchI];
            const tmp<scalarField>& Phii1 = psiE.boundaryField()[patchI].patchInternalField();
            Cb1 = Ci1/(1.0 + zPlus.value()*(Phib1 - Phii1));
        }
    }

    forAll(cMinus.boundaryField(),patchI)
    {
        if(cMinus.boundaryField()[patchI].type() == "fixedIonicFlux")
        {
            scalarField& Cb2 = cMinus.boundaryField()[patchI];
            const tmp<scalarField>& Ci2 = cMinus.boundaryField()[patchI].patchInternalField();
            scalarField& Phib2 = psiE.boundaryField()[patchI];
            const tmp<scalarField>&  Phii2 = psiE.boundaryField()[patchI].patchInternalField();
            Cb2 = Ci2/(1.0 + zMinus.value()*(Phib2 - Phii2));
        }
    }

    //#include "boundC.H"
    scalar residual;
    residual = Foam::sqrt(lengthcPlus + lengthcMinus + lengthpsiE);
    if (residual < convergenceCriterionPNP){ break; }
}

