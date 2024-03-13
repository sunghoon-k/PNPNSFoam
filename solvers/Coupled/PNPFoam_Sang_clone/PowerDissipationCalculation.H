word PowerDissipationCalculation(electrolyteProperties.lookup("PowerDissipationCalculation"));
if (PowerDissipationCalculation == "yes")
{
    scalar D1_ = D1.value();
    scalar D2_ = D2.value();
    scalar z1_ = z1.value();
    scalar z2_ = z2.value();
    scalar F_ = F.value();
    scalar T_ = T.value();
    scalar R_ = R.value();
    scalar C0_ = C0.value();
    scalar Phi0_ = Phi0.value();
    scalar U0_ = U0.value();
    scalar p0_ = p0.value();

    volScalarField P_total,P_electrical,P1,P2,P3;

    surfaceScalarField C1f = fvc::interpolate(C1*C0_);
    surfaceScalarField C2f = fvc::interpolate(C2*C0_);
    surfaceScalarField gradPhif = fvc::snGrad(Phi*Phi0_);
    surfaceScalarField magSf = mag(mesh.Sf());
    surfaceScalarField pf = fvc::interpolate(p*p0_);
    surfaceScalarField Phif = fvc::interpolate(Phi*Phi0_);
    volVectorField gradC1 = fvc::grad(C1*C0_);
    volVectorField gradC2 = fvc::grad(C2*C0_);
    volVectorField gradPhi = fvc::grad(Phi*Phi0_);

    // Ohmic Dissipation
    P1 = -(sqr(F_)/R_/T_)*(D1_*sqr(z1_)*C1*C0_ + D2_*sqr(z2_)*C2*C0_)*fvc::surfaceIntegrate(gradPhif*magSf)*cV;
    // Dissipation due to formation of extended space charge layer
    P2 = -F_*((z1_*D1_*gradC1 + z2_*D2_*gradC2) & gradPhi)*cV ;
    // Viscous dissipation
    dimensionedScalar dimp("dimp",dimensionSet(0,-2,2,0,0,0,0), 1);
    dimensionedScalar dimP3("dimP3",dimensionSet(0,-2,1,0,0,0,0), 1);
    P3 = F_*(z1_*C1*C0_ + z2_*C2*C0_)*fvc::surfaceIntegrate(U0_*phi*Phif)*cV*dimP3 + fvc::surfaceIntegrate(U0_*phi*pf*dimp)*cV*dimP3;
    // Electrical dissipation
    P_electrical = P1 + P2;
    // Total power dissipation
    P_total = P1 + P2 + P3;
}
