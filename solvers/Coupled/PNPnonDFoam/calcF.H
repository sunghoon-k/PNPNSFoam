volScalarField cV ( IOobject ( "cV", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ), mesh, dimensionedScalar("zero", dimensionSet(0, 3, 0, 0, 0, 0, 0), 0.0) );   

scalarField& cV_ = cV.internalField();        

forAll(cV_, k) {cV_[k] = mesh.V()[k];}        

surfaceScalarField jPlus =  -DPlus_nonD *(gradcPlusf  + zPlus*cPlusf*gradpsiEf  )*magSf; //+Pe*cPlusf*phi/(dimphi*l0_one))*magSf; // -DPlus*(gradcPlusf + zPlus*cPlusf*gradpsiEf/psiE0)*magSf; // [1/s]
surfaceScalarField jMinus = -DMinus_nonD*(gradcMinusf + zMinus*cMinusf*gradpsiEf)*magSf; // +Pe*cMinusf*phi/(dimphi*l0_one))*magSf; // -DMinus*(gradcMinusf + zMinus*cMinusf*gradpsiEf/psiE0)*magSf; // [1/s]        

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

scalar sqrDebL_nonD = 1/sqr(DebL_nonD.value());          
volScalarField zPluscPlusVol = zPlus*cPlus*(sqrDebL_nonD/2); // e*zPlus*cPlus/(eps0*epsr)*cV; // cPlus source of psiE's equation 
volScalarField zMinuscMinusVol = zMinus*cMinus*(sqrDebL_nonD/2); // e*zMinus*cMinus/(eps0*epsr)*cV; // cMinus source of psiE's equation

volScalarField zPluscPlusfVol = zPlus*fvc::surfaceIntegrate(cPlusf)*(sqrDebL_nonD/2)*cV;
volScalarField zMinuscMinusfVol = zMinus*fvc::surfaceIntegrate(cMinusf)*(sqrDebL_nonD/2)*cV;

volScalarField FcPlus = (fvc::laplacian(-DPlus_nonD, cPlus) + fvc::div(DPluszPlusphiE, cPlus) - fvc::surfaceIntegrate(jPlus))*cV;
volScalarField FcMinus = (fvc::laplacian(-DMinus_nonD, cMinus) + fvc::div(DMinuszMinusphiE, cMinus) - fvc::surfaceIntegrate(jMinus))*cV*sqr(l0_one);       
volScalarField FpsiE = (fvc::laplacian(psiE)*sqr(l0_one) + zPluscPlusVol + zMinuscMinusVol - fvc::surfaceIntegrate(gradpsiEf*magSf)*sqr(l0_one) -zPluscPlusfVol -zMinuscMinusfVol)*cV;            

// Assign F to blockM's right-hand side            
forAll(PNPEqn.source(),cellI)            
{                
    PNPEqn.source()[cellI](0) += FcPlus[cellI];                
    PNPEqn.source()[cellI](1) += FcMinus[cellI];                
    PNPEqn.source()[cellI](2) += FpsiE[cellI];            
}            
