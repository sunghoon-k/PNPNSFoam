runTime.setTime(runTime.startTime(), 0);
scalar Isum(0);
scalar imin,imax,iinst,istep(0);
scalar ipnp_ns2,dC2(10);
while(runTime.loop())
{
    #include "readBlockSolverControls.H"
    ipnp_ns2 = 0;
    dC2=dCmin+10;
    while (ipnp_ns2<nNSPNPIteration and dC2>dCmin)
    {
        volScalarField C1old = C1;
        Info<<dVstring<<phiInstant_<<"           PNP-NS/time "<<runTime.timeName()<<"/iteration #"<<ipnp_ns2<<"/<Transient Poisson-Nernst-Planck>"<<endl<<endl;
        #include "PNPEqn_t.H"
        Info<<endl;
        Info<<dVstring<<phiInstant_<<"           PNP-NS/time "<<runTime.timeName()<<"/iteration #"<<ipnp_ns2<<"/<Navier-Stokes>"<<endl<<endl;
        #include "NSEqn.H"
        Info<<endl;
        dC2 = mag(gMax((C1.internalField() - C1old.internalField())));
        Info <<nl<<dVstring<<phiInstant_<<": C convergence: "<<dC2<<endl<<endl;
        ipnp_ns2++;
    }
    desalt.CurrentPerTimestep(runTime.timeName(),istep==0);
    if (istep > 5)
    Isum += desalt.CurrentSum();
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;
    runTime.write();
    istep++;
}
Info <<endl;
runTime.write();

// ************************************************************************* //
