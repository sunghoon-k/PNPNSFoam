#ifndef PNP_NSEQN_H_INCLUDED
#define PNP_NSEQN_H_INCLUDED


scalar nNSPNPIteration(readScalar(runTime.controlDict().lookup("nNSPNPIteration")));
scalar Isum(0);
scalarList Is(nNSPNPIteration,0);
scalar ipnp_ns,dC(10);
scalar dCmin(readScalar(runTime.controlDict().lookup("C1C2convergence")));
word dVstring(" >>> dV = ");
{
  #include "readBlockSolverControls.H"
  ipnp_ns = 0;dC=dCmin+10;
  while(ipnp_ns < nNSPNPIteration and dC>dCmin)
  {
    if(!solvedTransient)
    {
      Info <<dVstring<<phiInstant_<<":         Steadystate solver/PNP-NS Iteration #" << ipnp_ns <<"/<Steadystate Poisson-Nernst-Planck>  "<<endl<<endl;
      ///C1.oldTime();
      volScalarField C1old = C1;
      #include "PNPEqn.H"
      /// Solve NS Equations:
      Info <<endl<<endl;
      Info <<dVstring<<phiInstant_<<":         Steadystate solver/PNP-NS Iteration #" << ipnp_ns <<"/<Steadystate Navier-Stokes>  "<<endl<<endl;
      #include "NSEqn.H"
      #include "NSEqn.H"
      #include "NSEqn.H"
      Info <<endl;
      dC = mag(gMax((C1.internalField() - C1old.internalField())));
      Info <<nl<<dVstring<<phiInstant_<<": C convergence: "<<dC<<endl<<endl;
      ipnp_ns++;
    }
    if((ipnp_ns==nNSPNPIteration and dC>dCmin) or solvedTransient)
    {
      Info <<dVstring<<phiInstant_<<":         !!!!! Starting transient solver  "<<endl<<endl;
      #include "transient_PNP-NSDImlessFoam.H"
      solvedTransient = true;
      break;
    }
  }
}

#endif
