/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    blockCoupledScalarTransportFoam

Description
    Solves two coupled transport equations in a block-coupled manner

        1) transport equation for a passive scalar
        2) diffusion only

    This resembles heat exchanging flow through a porous medium

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvBlockMatrix.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "ElectrochemicalSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "NormalizedParameters.H"

    word PNPsolverType(mesh.solutionDict().subDict("solvers").subDict("blockPNP").lookup("solver"));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport  \n" << endl;
    bool solveTransient(readBool(runTime.controlDict().lookup("PNPNStransient")));
    bool solvedTransientInSteadyStateMode = false;
    bool continueSolve = false;
    double timeEnd = runTime.endTime().value();
    double timeLap = 1.0;
    scalar phiInstant_ = ECsystem.phiInstant();


    forAll(mesh.boundary(),patchI)
    {
        if(mesh.boundary()[patchI].name() == phiDyMBound)
        {
            scalarField& Phib = Phi.boundaryField()[patchI];
            if(Phib[0]-phiInstant_>0.01)
              { ECsystem.set_phiInstant(Phib[0]); continueSolve = true; break;}
        }
    }


    while(ECsystem.phiRun())
    {
        phiInstant_ = ECsystem.phiInstant();
        Info << "\n-------------- | Phi = " << phiInstant_<< " |--------------"<<endl;
        forAll(mesh.boundary(),patchI)
        {
            if(mesh.boundary()[patchI].name() == phiDyMBound)
            {
                scalarField& Phib = Phi.boundaryField()[patchI];
                forAll(Phib, faceI)
                {
                    Phib[faceI] = phiInstant_;
                }
            }
        }


      runTime.setEndTime(timeEnd*timeLap);
      scalar nNSPNPIteration(readScalar(runTime.controlDict().lookup("nNSPNPIteration")));
      scalarList Is(nNSPNPIteration,0);
      scalar ipnp_ns = 0,dC(10);
      scalar dCmin(readScalar(runTime.controlDict().lookup("C1C2convergence")));
      word dVstring(" >>> dV = ");
      char phiInstant_str[6];
      gcvt(phiInstant_,6,phiInstant_str);

      while(runTime.loop())
      {
        #include "readBlockSolverControls.H"
        ipnp_ns = 0;dC=dCmin+10;

        if(solveTransient or solvedTransientInSteadyStateMode)
        {
          Info <<dVstring<<phiInstant_<<":         Transient solver/PNP-NS Step      # " << runTime.timeIndex() <<endl;
          C1.storePrevIter();
          C2.storePrevIter();
          Phi.storePrevIter();
          while(ipnp_ns < nNSPNPIteration and dC>dCmin)
          {
            Info <<"         Transient solver/PNP-NS Iteration      # " << ipnp_ns <<endl;
            ///C1.oldTime();
            volScalarField C1old = C1;
            Info<<endl<<"*** PNPEqn_t ***"<<endl;
            #include "PNPEqn_t.H"
            dC = mag(gMax((C1.internalField() - C1old.internalField())));
            Info <<endl<<nl<<dVstring<<phiInstant_<<": C convergence: "<<dC<<endl<<endl;
            ipnp_ns++;
          }
          if( transportProperties.subDict("Desalination").lookupOrDefault<Switch>("saveCurrentPerTimestep",false))
            ECsystem.CurrentPerTimestep(runTime.timeIndex());
        }
        else
        {
          Info <<dVstring<<phiInstant_<<":         SteadyState solver"<<endl;
          while(ipnp_ns < nNSPNPIteration and dC>dCmin)
          {
            Info <<"         SteadyState solver/PNP-NS Iteration      # " << ipnp_ns <<endl;
            volScalarField C1old = C1;
            Info<<endl<<"*** PNPEqn ***"<<endl;
            #include "PNPEqn.H"
            dC = mag(gMax((C1.internalField() - C1old.internalField())));
            Info <<endl<<nl<<dVstring<<phiInstant_<<": C convergence: "<<dC<<endl<<endl;
            ipnp_ns++;
            if(ipnp_ns==nNSPNPIteration and dC>dCmin)
            {
              solvedTransientInSteadyStateMode = true; /// after this, the solver will run into the transient mode to finish this applied Voltage
            }
          }
        }
        runTime.write();
        if (!solveTransient and !solvedTransientInSteadyStateMode) break;
      }

      timeLap+=1.0;
      ECsystem.changePhi();

    }

    return 0;
}






// ************************************************************************* //
