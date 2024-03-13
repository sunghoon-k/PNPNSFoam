Info << "****************************************************************\n" << endl;
Info << "Scale Parameters" << endl;
Info << "****************************************************************" << endl;
//Calculate reference values for DimensionLess-Form of Equations
// Calculate average-diffusivity

dimensionedScalar l0(transportProperties.subDict("Electrolyte").lookup("l0"));
dimensionedScalar C0(transportProperties.subDict("Electrolyte").lookup("C0"));
dimensionedScalar T(transportProperties.subDict("Electrolyte").lookup("T"));
dimensionedScalar e(transportProperties.subDict("Electrolyte").lookup("e"));
dimensionedScalar kB(transportProperties.subDict("Electrolyte").lookup("kB"));
dimensionedScalar F(transportProperties.subDict("Electrolyte").lookup("F"));
dimensionedScalar NA(transportProperties.subDict("Electrolyte").lookup("NA"));
dimensionedScalar eps0(transportProperties.subDict("Electrolyte").lookup("eps0"));
dimensionedScalar epsr(transportProperties.subDict("Electrolyte").lookup("epsr"));
// dimensionedScalar wsc(transportProperties.subDict("Electrolyte").lookup("wsc"));
dimensionedScalar rho(transportProperties.subDict("Electrolyte").lookup("rho"));
dimensionedScalar mu(transportProperties.subDict("Electrolyte").lookup("mu"));
dimensionedScalar nu(transportProperties.subDict("Electrolyte").lookup("nu"));
//dimensionedScalar mu = nu/rho;
dimensionedScalar z1(transportProperties.subDict("Ions").subDict("C1").lookup("z"));
dimensionedScalar D1(transportProperties.subDict("Ions").subDict("C1").lookup("D"));
dimensionedScalar z2(transportProperties.subDict("Ions").subDict("C2").lookup("z"));
dimensionedScalar D2(transportProperties.subDict("Ions").subDict("C2").lookup("D"));
// dimensionedScalar D0("D0",dimensionSet(0,2,-1,0,0,0,0), (D1.value() + D2.value())/2);
dimensionedScalar D0("D0", D1 + D2);
Info << "Average-diffusivity: " << D0.value() << endl;

 // Calculate time-scale
dimensionedScalar t0("t0",dimensionSet(0,0,1,0,0,0,0), sqr(l0.value())/D0.value());
Info << "Time-scale: " << t0.value() << endl;

 // Calculate potential-scale(Thermal voltage)
dimensionedScalar Phi0("Phi0",dimensionSet(1,2,-3,0,0,-1,0), (kB.value()*T.value())/(z1.value()*e.value()));
Info << "Potential-scale: " << Phi0.value() << endl;

 // Calculate IonicFluxDensity-scale
// dimensionedScalar IFD0("IFD0",dimensionSet(0,-3,0,0,1,0,0), D0.value()*C0.value()/l0.value());
// Info << "IonicFluxDensity-scale: " << IFD0.value() << endl;

 // Calculate velocity-scale
dimensionedScalar U0("U0",dimensionSet(0,1,-1,0,0,0,0), epsr.value()*eps0.value()*sqr(Phi0.value())/(mu.value()*l0.value()));
Info << "Velocity-scale: " << U0.value() << endl;

 // Calculate pressure-scale
dimensionedScalar p0("p0",dimensionSet(1,-1,-2,0,0,0,0), mu.value()*U0.value()/l0.value());
Info << "pressure-scale: " << p0.value() << endl;

 // Calculate surfacecharge-scale
dimensionedScalar SC("SC",dimensionSet(0,-2,1,0,0,1,0), e.value()*C0.value()*l0.value()*NA.value());
Info << "surfacecharge-scale: " << SC.value() << endl;

 // Calculate Debye-Length
dimensionedScalar DebL("DebL",dimensionSet(0,1,0,0,0,0,0), Foam::sqrt(epsr.value()*eps0.value()*kB.value()*T.value()/C0.value()/NA.value()/sqr(e.value())/2));
Info << "Debye-Length: " << DebL.value() << endl;
 // Calculate dimensionless diffusion coefficients
dimensionedScalar DimlD1("DimlD1",dimensionSet(0,0,0,0,0,0,0), D1.value()/D0.value());

Info << "Dimensionless diffusion coefficient DimlessD1: " << DimlD1.value() << endl;

dimensionedScalar DimlD2("DimlD2",dimensionSet(0,0,0,0,0,0,0), D2.value()/D0.value());

Info << "Dimensionless diffusion coefficient DimlessD2: " << DimlD2.value() << endl;

 // Calculate dimensionless Debye-Length
scalar Epsilon = DebL.value()/l0.value();
scalar Epsilon2 = sqr(Epsilon);
Info << "Dimensionless Debye-Length: " << Epsilon << endl;
Info << "Epsilon2: " << Epsilon2 << endl;

 // Calculate Schmidt Number
scalar Sc = mu.value()/(rho.value()*D0.value()) ;
Info << "Schmidt Number: " << Sc << endl;

 // Calculate Reynold Number
scalar Re = U0.value()*l0.value()*rho.value()/mu.value() ;
Info << "Reynold Number: " << Re << endl;

 // Calculate Pelec Number
scalar Pe = U0.value()*l0.value()/D0.value() ;
Info << "Pelec Number: " << Pe << endl;

Info << "****************************************************************" << endl;
