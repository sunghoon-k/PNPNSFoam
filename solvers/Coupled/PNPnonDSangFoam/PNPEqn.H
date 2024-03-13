reachedResidual = false;

/********************************************************************************************************/
volScalarField cPlusT(cPlus);        
volScalarField cMinusT(cMinus);        
volScalarField psiET(psiE);        
volScalarField cV ( IOobject ( "cV", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ), mesh, dimensionedScalar("zero", dimensionSet(0, 3, 0, 0, 0, 0, 0), 0.0) );   

scalarField& cV_ = cV.internalField();        

forAll(cV_, k) {cV_[k] = mesh.V()[k];}        
/********************************************************************************************************/

scalar nNewtonIter(readScalar(runTime.controlDict().lookup("nNewtonIter")));
scalar convergenceCriterion(readScalar(mesh.solutionDict().subDict("solvers").subDict("PNP").lookup("convergenceCriterion")));

for(int NewtonIter = 0; NewtonIter < nNewtonIter; NewtonIter++)
{

    surfaceScalarField gradcPlusf   = fvc::snGrad(cPlus);            
    surfaceScalarField gradcMinusf  = fvc::snGrad(cMinus);            
    surfaceScalarField cPlusf       = fvc::interpolate(cPlus);            
    surfaceScalarField cMinusf      = fvc::interpolate(cMinus);            
    surfaceScalarField gradpsiEf    = fvc::snGrad(psiE);            
    surfaceScalarField magSf        = mag(mesh.Sf());            

    // Calculate J-Flux            
    surfaceScalarField jPlus =  -DPlus_nonD *(gradcPlusf  + zPlus*cPlusf*gradpsiEf)*magSf; // -DPlus*(gradcPlusf + zPlus*cPlusf*gradpsiEf/psiE0)*magSf; // [1/s]
    surfaceScalarField jMinus = -DMinus_nonD*(gradcMinusf + zMinus*cMinusf*gradpsiEf)*magSf; // -DMinus*(gradcMinusf + zMinus*cMinusf*gradpsiEf/psiE0)*magSf; // [1/s]        

    // Assign J-Flux to specific value            
    forAll(cPlus.boundaryField(),patchI)            
    {                
    if(cPlus.boundaryField()[patchI].type() == "fixedIonicFlux")                
    {                    
        scalarField& jPlusf = jPlus.boundaryField()[patchI];                    
        forAll(jPlusf, i)                    
        {                        
            jPlusf[i] = 0;                    
        }                
    }            
    }            

    // Assign J-Flux to specific value            
    forAll(cMinus.boundaryField(),patchI)            
    {                
    if(cMinus.boundaryField()[patchI].type() == "fixedIonicFlux")                
    {                    
        scalarField& jMinusf = jMinus.boundaryField()[patchI];                    
        forAll(jMinusf, i)                    
        {                        
            jMinusf[i] = 0;                    
        }                
    }            
    }            

    // Calculate Right-Hand Side F  
    scalar sqrDebL_nonD = 1/sqr(DebL_nonD.value());          
    volScalarField zPluscPlusVol = zPlus*cPlus*(sqrDebL_nonD/2)*cV; // e*zPlus*cPlus/(eps0*epsr)*cV; // cPlus source of psiE's equation 
    volScalarField zMinuscMinusVol = zMinus*cMinus*(sqrDebL_nonD/2)*cV; // e*zMinus*cMinus/(eps0*epsr)*cV; // cMinus source of psiE's equation

    volScalarField FcPlus = fvc::surfaceIntegrate(jPlus)*cV*sqr(l0_one);
    volScalarField FcMinus = fvc::surfaceIntegrate(jMinus)*cV*sqr(l0_one);       
    volScalarField FpsiE = fvc::surfaceIntegrate(gradpsiEf*magSf)*cV*sqr(l0_one) + (zPluscPlusVol + zMinuscMinusVol);            

    bool bounded = false;

    surfaceScalarField DPluszPlusphiE("DPluszPlusphiE",             -DPlus_nonD*zPlus*gradpsiEf*magSf); // 
    surfaceScalarField DMinuszMinusphiE("DMinuszMinusphiE",         -DMinus_nonD*zMinus*gradpsiEf*magSf);
    surfaceScalarField DPluszPluscPlusf("DPluszPluscPlusf",         -DPlus_nonD*zPlus*cPlusf);
    surfaceScalarField DMinuszMinuscMinusf("DMinuszMinuscMinusf",   -DMinus_nonD*zMinus*cMinusf);
    
    forAll(cPlus.boundaryField(),patchI)
    {
        if(cPlus.boundaryField()[patchI].type() == "fixedIonicFlux")
        {
            scalarField& DPluszPlusphiE_ = DPluszPlusphiE.boundaryField()[patchI];
            forAll(DPluszPlusphiE_, i) { DPluszPlusphiE_[i] = 0; }
            scalarField& DPluszPluscPlusf_ = DPluszPluscPlusf.boundaryField()[patchI];
            forAll(DPluszPluscPlusf_, i) { DPluszPluscPlusf_[i] = 0; }
        }
    }

    forAll(cMinus.boundaryField(),patchI)
    {
        if(cMinus.boundaryField()[patchI].type() == "fixedIonicFlux")
        {
            //Info<<nl<<"********************************zero!!!"<<endl;

            scalarField& DMinuszMinusphiE_ = DMinuszMinusphiE.boundaryField()[patchI];
            forAll(DMinuszMinusphiE_, i) { DMinuszMinusphiE_[i] = 0; }
            scalarField& DMinuszMinuscMinusf_ = DMinuszMinuscMinusf.boundaryField()[patchI];
            forAll(DMinuszMinuscMinusf_, i) { DMinuszMinuscMinusf_[i] = 0; }
        }
    }



    /******************************** couplingTerms.H ********************************/
    fvBlockMatrix<vector3> PNPEqn(PNP);
    #include "cEqn.H"
    #include "psiEEqn.H"
    #include "couplingTerms.H"

    /****************************************************************/

    // Assign F to blockM's right-hand side            
    forAll(PNPEqn.source(),cellI)            
    {                
        PNPEqn.source()[cellI](0) = -FcPlus[cellI];                
        PNPEqn.source()[cellI](1) = -FcMinus[cellI];                
        PNPEqn.source()[cellI](2) = -FpsiE[cellI];            
    }            

    maxResidual = cmptMax(PNPEqn.solve().initialResidual()); // residual = b - A*x_n
    
    // Retrieve solution
    PNPEqn.retrieveSolution(0, cPlus.internalField());
    PNPEqn.retrieveSolution(1, cMinus.internalField());
    PNPEqn.retrieveSolution(2, psiE.internalField());

    scalar lengthcPlus(0);
    scalar lengthcMinus(0);
    scalar lengthpsiE(0);

    // cPlus residual-length

    forAll(cPlus.internalField(),i)
    {
        lengthcPlus += sqr(cPlus.internalField()[i]);
    }

    // cMinus residual-length

    forAll(cMinus.internalField(),i)
    {
        lengthcMinus += sqr(cMinus.internalField()[i]);
    }

    // Phi residual-length

    forAll(psiE.internalField(),i)
    {
        lengthpsiE += sqr(psiE.internalField()[i]);
    }

    Info <<" Iter # "<< PNPIter  << " - Residual: " <<"[" << "cPlus:" << Foam::sqrt(lengthcPlus)
    << "  " << "cMinus:" << Foam::sqrt(lengthcMinus) << "  " << "Phi:" << Foam::sqrt(lengthpsiE) << "]"  << endl;

    // Correct cPlus,cMinus,Phi

    cPlusT  += cPlus;
    cMinusT  += cMinus;
    psiET += psiE;
    // Assign cPlusT,cMinusT,PhiT to cPlus,cMinus,Phi
    cPlus  = cPlusT;
    cMinus  = cMinusT;
    psiE = psiET;

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

    #include "boundC.H"

    scalar residual;
    residual = Foam::sqrt(lengthcPlus + lengthcMinus + lengthpsiE);
    if (residual < convergenceCriterion)
    {
        break;
    }
    //#include "convergenceCheck.H"
    netCharge = (zPlus*cPlus + zMinus*cMinus);
}

//Info<< "maxResidual = " << maxResidual << nl
//    << endl;
