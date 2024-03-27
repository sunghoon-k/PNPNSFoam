surfaceScalarField presSource( "presSource", fvc::interpolate(rAU)*((fvc::interpolate(fvc::grad(p)))& mesh.Sf()));

volVectorField FU = (Re * fvc::div(phi, U) - fvc::laplacian(nu_one, U) + fvc::grad(p))*cV;
volScalarField Fp = (fvc::div(phi) - fvc::laplacian(rAU, p) + fvc::div(presSource))*cV;
// Assign F to blockM's right-hand side            
forAll(UpEqn.source(),cellI)            
{                
    UpEqn.source()[cellI](0) = -FU[cellI](0);                
    UpEqn.source()[cellI](1) = -FU[cellI](1);                
    UpEqn.source()[cellI](2) = -FU[cellI](2);                
    UpEqn.source()[cellI](3) = -Fp[cellI];                
}            
