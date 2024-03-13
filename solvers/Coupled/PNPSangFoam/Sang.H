surfaceScalarField gradcPlusf = fvc::snGrad(cPlus);            
surfaceScalarField gradcMinusf = fvc::snGrad(cMinus);            
surfaceScalarField cPlusf = fvc::interpolate(cPlus);            
surfaceScalarField cMinusf = fvc::interpolate(cMinus);            
surfaceScalarField gradpsiEf = fvc::snGrad(psiE);            
surfaceScalarField magSf = mag(mesh.Sf());            

// Calculate J-Flux            
surfaceScalarField jPlus = -DPlus*(gradcPlusf + zPlus*cPlusf*gradpsiEf/psiE0)*magSf; // [1/s]
surfaceScalarField jMinus = -DMinus*(gradcMinusf + zMinus*cMinusf*gradpsiEf/psiE0)*magSf; // [1/s]        

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
volScalarField zPluscPlusVol = e*zPlus*cPlus/(eps0*epsr)*cV; // cPlus source of psiE's equation 
volScalarField zMinuscMinusVol = e*zMinus*cMinus/(eps0*epsr)*cV; // cMinus source of psiE's equation

// Info << jPlus.dimensions() << endl; // dimension = D * c * l0 = 1/s
volScalarField FcPlus = fvc::surfaceIntegrate(jPlus)*cV;
// Info << (fvc::surfaceIntegrate(jPlus)).dimensions() << endl; // dimension = D * c / l0^2 = 1/(m^3 * s)
// Info << FcPlus.dimensions() << endl; // dimension = D * c * l0 = 1/s

volScalarField FcMinus = fvc::surfaceIntegrate(jMinus)*cV;       
volScalarField FpsiE = fvc::surfaceIntegrate(gradpsiEf*magSf)*cV + (zPluscPlusVol + zMinuscMinusVol);            
