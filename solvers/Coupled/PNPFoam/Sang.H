surfaceScalarField gradcPlusf = fvc::snGrad(cPlus);            
surfaceScalarField gradcMinusf = fvc::snGrad(cMinus);            
surfaceScalarField cPlusf = fvc::interpolate(cPlus);            
surfaceScalarField cMinusf = fvc::interpolate(cMinus);            
surfaceScalarField gradpsiEf = fvc::snGrad(psiE);            
surfaceScalarField magSf = mag(mesh.Sf());            
