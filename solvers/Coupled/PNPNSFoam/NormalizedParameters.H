Info << "****************************************************************\n" << endl;
Info << "Scale Parameters" << endl;
Info << "****************************************************************" << endl;
//Calculate reference values for DimensionLess-Form of Equations
// Calculate average-diffusivity

dimensionedScalar l0(ECsystem.l0);
dimensionedScalar c0(ECsystem.c0);
dimensionedScalar T(ECsystem.T);
dimensionedScalar e(ElectrochemicalSystem::e);
dimensionedScalar kB(ElectrochemicalSystem::kB);
dimensionedScalar F(ElectrochemicalSystem::F);
dimensionedScalar NA(ElectrochemicalSystem::NA);
dimensionedScalar eps0(ElectrochemicalSystem::eps0);
dimensionedScalar epsr(ECsystem.epsr);
// dimensionedScalar wsc(transportProperties.subDict("Electrolyte").lookup("wsc"));
dimensionedScalar rho(ECsystem.rho);
// dimensionedScalar mu(transportProperties.subDict("Electrolyte").lookup("mu"));
dimensionedScalar mu(ECsystem.mu);
dimensionedScalar nu = mu/rho;
dimensionedScalar zPlus(ECsystem.zPlus);
dimensionedScalar DPlus(ECsystem.DPlus);
dimensionedScalar zMinus(ECsystem.zMinus);
dimensionedScalar DMinus(ECsystem.DMinus);
// dimensionedScalar D0("D0",dimensionSet(0,2,-1,0,0,0,0), (D1.value() + D2.value())/2);
dimensionedScalar D0("D0", (DPlus + DMinus)/2);
Info << "Average-diffusivity: " << D0.value() << endl;

 // Calculate time-scale
// dimensionedScalar t0("t0",dimensionSet(0,0,1,0,0,0,0), sqr(l0.value())/D0.value());
dimensionedScalar t0("t0", sqr(l0)/D0);
Info << "Time-scale: " << t0.value() << endl;

 // Calculate potential-scale(Thermal voltage)
// dimensionedScalar psiE0("psiE0",dimensionSet(1,2,-3,0,0,-1,0), (kB_const.value()*T_const.value())/(zPlus.value()*e_const.value()));
dimensionedScalar psiE0("psiE0", (kB*T)/(zPlus*e));
Info << "Potential-scale: " << psiE0.value() << endl;

 // Calculate IonicFluxDensity-scale
// dimensionedScalar IFD0("IFD0",dimensionSet(0,-3,0,0,1,0,0), D0.value()*C0.value()/l0.value());
// Info << "IonicFluxDensity-scale: " << IFD0.value() << endl;

 // Calculate velocity-scale
dimensionedScalar U0("U0", epsr*eps0*sqr(psiE0)/(nu*l0));
Info << "Velocity-scale: " << U0.value() << endl;

 // Calculate pressure-scale
dimensionedScalar p0("p0", nu*U0/l0);
Info << "pressure-scale: " << p0.value() << endl;

 // Calculate surfacecharge-scale
dimensionedScalar SC("SC", e*c0*l0*NA);
Info << "surfacecharge-scale: " << SC.value() << endl;

 // Calculate Debye-Length
dimensionedScalar DebL("DebL", Foam::sqrt(epsr*eps0*kB*T/(2*c0*sqr(e))));
Info << "Debye-Length: " << DebL.value() << endl;

 // Calculate dimensionless diffusion coefficients
dimensionedScalar DimlDPlus("DimlDPlus", DPlus/D0);

Info << "Dimensionless diffusion coefficient DimlessDPlus: " << DimlDPlus.value() << endl;

dimensionedScalar DimlDMinus("DimlDMinus", DMinus/D0);

Info << "Dimensionless diffusion coefficient DimlessDMinus: " << DimlDMinus.value() << endl;

 // Calculate dimensionless Debye-Length
scalar Epsilon = (DebL/l0).value();
scalar Epsilon2 = sqr(Epsilon);
Info << "Dimensionless Debye-Length: " << Epsilon << endl;
Info << "Epsilon2: " << Epsilon2 << endl;

 // Calculate Schmidt Number
scalar Sc = (nu/D0).value();
Info << "Schmidt Number: " << Sc << endl;

 // Calculate Reynold Number
scalar Re = (U0*l0/nu).value();
Info << "Reynold Number: " << Re << endl;

 // Calculate Pelec Number
scalar Pe = (U0*l0/D0).value();
Info << "Pelec Number: " << Pe << endl;

Info << "****************************************************************" << endl;

/*

*/
