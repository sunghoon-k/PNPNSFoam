Info << "****************************************************************\n" << endl;
Info << "Scale Parameters" << endl;
Info << "****************************************************************" << endl;
//Calculate reference values for DimensionLess-Form of Equations
// Calculate average-diffusivity

dimensionedScalar l0(ECsystem.l0);
dimensionedScalar l0_one("l0_one", l0.dimensions(), 1.0);
dimensionedScalar dimphi("dimphi",dimensionSet(0,3,-1,0,0,0,0), 1);
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
dimensionedScalar D0("D0", DPlus + DMinus);
//Info << "Average-diffusivity: " << D0.value() << endl;
Info << "Average-diffusivity's value: " << D0.value() << " " << D0.dimensions()<< endl;

 // Calculate time-scale
// dimensionedScalar t0("t0",dimensionSet(0,0,1,0,0,0,0), sqr(l0.value())/D0.value());
dimensionedScalar t0("t0", sqr(l0)/D0);
// Info << "Time-scale: " << t0.value() << endl;
Info << "Time-scale's value: " << t0.value() << " " << t0.dimensions()<< endl;
dimensionedScalar t0_one("t0_one", t0.dimensions(), 1.0);

 // Calculate potential-scale(Thermal voltage)
// dimensionedScalar psiE0("psiE0",dimensionSet(1,2,-3,0,0,-1,0), (kB_const.value()*T_const.value())/(zPlus.value()*e_const.value()));
dimensionedScalar psiE0("psiE0", (kB*T)/(e));
// Info << "Potential-scale: " << psiE0.value() << endl;
Info << "Potential-scale's value: " << psiE0.value() << " " << psiE0.dimensions()<< endl;

 // Calculate IonicFluxDensity-scale
// dimensionedScalar IFD0("IFD0",dimensionSet(0,-3,0,0,1,0,0), D0.value()*C0.value()/l0.value());
// Info << "IonicFluxDensity-scale: " << IFD0.value() << endl;

 // Calculate velocity-scale
dimensionedScalar U0("U0", epsr*eps0*sqr(psiE0)/(mu*l0));
// Info << "Velocity-scale: " << U0.value() << endl;
Info << "Velocity-scale's value: " << U0.value() << " " << U0.dimensions()<< endl;

 // Calculate pressure-scale
dimensionedScalar p0("p0", nu*U0/l0);
// Info << "pressure-scale: " << p0.value() << endl;
Info << "pressure-scale's value: " << p0.value() << " " << p0.dimensions()<< endl;

 // Calculate surfacecharge-scale
// dimensionedScalar SC("SC", e*c0*l0*NA);
// Info << "surfacecharge-scale: " << SC.value() << endl;

 // Calculate Debye-Length
dimensionedScalar DebL("DebL", Foam::sqrt(epsr*eps0*kB*T/(2*c0*NA*sqr(e))));
Info << "Debye-Length's value: " << DebL.value() << " " << DebL.dimensions()<< endl;

 // Calculate dimensionless diffusion coefficients
dimensionedScalar DPlus_nonD("DPlus_nonD", (DPlus/D0).value());
// scalar DPlus_nonD = (DPlus/D0).value();
Info << "Dimensionless diffusion coefficient DPlus_nonD: " << DPlus_nonD.value() << endl;

dimensionedScalar DMinus_nonD("DMinus_nonD", (DMinus/D0).value());
// scalar DMinus_nonD = (DMinus/D0).value();
Info << "Dimensionless diffusion coefficient DMinus_nonD: " << DMinus_nonD.value() << endl;

 // Calculate dimensionless Debye-Length
//scalar DebL_nonD = (DebL/l0).value();
dimensionedScalar DebL_nonD("DebL_nonD", DebL/l0);

// scalar Epsilon2 = sqr(Epsilon);
Info << "Dimensionless Debye-Length: " << DebL_nonD.value() << endl;
// Info << "Epsilon2: " << Epsilon2 << endl;

scalar oneOver2SquareDebL_nonD = 1/(2*sqr(DebL_nonD.value()));

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
dimensionedScalar dimL("dimL",dimensionSet(0,1,0,0,0,0,0), 1);
dimensionedScalar dimU("dimU",dimensionSet(0,1,-1,0,0,0,0), 1);
dimensionedScalar nu_one("nu_one", nu.dimensions(), 1.0);
dimensionedScalar bodyForce_one("bodyForce_one", dimensionSet(0,2,-2,0,0,0,0), 1.0);
dimensionedScalar ddtc_one("ddtc_one", dimensionSet(0, -2, 1, 0, 0, 0, 0), 1.0);

/*

*/
