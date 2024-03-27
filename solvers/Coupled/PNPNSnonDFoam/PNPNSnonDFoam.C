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
    pUCoupledFoam

Description
    Steady-state solver for incompressible, turbulent flow, with implicit
    coupling between pressure and velocity achieved by fvBlockMatrix.
    Turbulence is in this version solved using the existing turbulence
    structure.

Authors
    Klas Jareteg, Chalmers University of Technology,
    Vuko Vukcevic, FMENA Zagreb.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvBlockMatrix.H"
#include "singlePhaseTransportModel.H"
//#include "RASModel.H"
//#include "simpleControl.H"
#include "ElectrochemicalSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    bool solveTransient(readBool(runTime.controlDict().lookup("PNPNStransient")));
    bool solveNS(readBool(runTime.controlDict().lookup("solveNS")));
    bool splitpsi(readBool(runTime.controlDict().lookup("splitpsi")));

//    simpleControl simple(mesh);
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createControls.H"
    #include "initConvergenceCheck.H"

    bool reachedResidual = false;
    Info<< "\nStarting time loop\n" << endl;

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        #include "readBlockSolverControls.H"
        #include "readFieldBounds.H"

        maxResidual = 10;

        reachedResidual = false;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if(solveTransient)
        {
            if(solveNS)
            {
                if(splitpsi){Info <<"         Transient solver/PNPNS_split Iteration      # " <<endl;}
                else{Info <<"         Transient solver/PNPNS Iteration      # " <<endl;}
            }
            else
            {
                if(splitpsi){Info <<"         Transient solver/PNP_split Iteration      # " <<endl;}
                else{Info <<"         Transient solver/PNP Iteration      # " <<endl;}
            }
        }
        else
        {
            if(solveNS)
            {
                if(splitpsi){Info <<"         Steady solver/PNPNS_split Iteration      # " <<endl;}
                else{Info <<"         Steady solver/PNPNS Iteration      # " <<endl;}
            }
            else
            {
                if(splitpsi){Info <<"         Steady solver/PNP_split Iteration      # " <<endl;}
                else{Info <<"         Steady solver/PNP Iteration      # " <<endl;}
            }
        }

        // scalar PNPIter = 0;
        // while(PNPIter++ < nPNPIter and maxResidual > 1e-6) // 
        for(int PNPIter = 0; PNPIter < nPNPIter; PNPIter++)
        {
            surfaceScalarField gradcPlusf   = fvc::snGrad(cPlus);            
            surfaceScalarField gradcMinusf  = fvc::snGrad(cMinus);            
            surfaceScalarField cPlusf       = fvc::interpolate(cPlus);            
            surfaceScalarField cMinusf      = fvc::interpolate(cMinus);            
            surfaceScalarField gradpsiEf    = fvc::snGrad(psiE);    
            surfaceScalarField gradpsiIf    = fvc::snGrad(psiI);    
            //surfaceVectorField intergradpsiEf = fvc::interpolate(fvc::grad(psiE));            
            surfaceScalarField magSf        = mag(mesh.Sf());            

            bool bounded = false;


            surfaceScalarField DPluszPlusphiE("DPluszPlusphiE",             -DPlus_nonD*zPlus*(gradpsiEf+gradpsiIf   + Pe*cPlusf*phi/(dimphi*l0_one))*magSf); // 
            surfaceScalarField DMinuszMinusphiE("DMinuszMinusphiE",         -DMinus_nonD*zMinus*(gradpsiEf+gradpsiIf + Pe*cMinusf*phi/(dimphi*l0_one))*magSf);
            surfaceScalarField DPluszPluscPlusf("DPluszPluscPlusf",         -DPlus_nonD*zPlus*cPlusf);
            surfaceScalarField DMinuszMinuscMinusf("DMinuszMinuscMinusf",   -DMinus_nonD*zMinus*cMinusf);
            
            forAll(cPlus.boundaryField(),patchI)
            {
                if(cPlus.boundaryField()[patchI].type() == "zeroIonicFlux_nonD")//  or "zeroIonicFlux_nonD") // "zeroIonicFlux_nonD" or
                {
                    //scalarField& DPluszPlusinterphiE_ = DPluszPlusinterphiE.boundaryField()[patchI];
                    //forAll(DPluszPlusinterphiE_, i) { DPluszPlusinterphiE_[i] = 0; }
                    scalarField& DPluszPlusphiE_ = DPluszPlusphiE.boundaryField()[patchI];
                    forAll(DPluszPlusphiE_, i) { DPluszPlusphiE_[i] = 0; }
                    scalarField& DPluszPluscPlusf_ = DPluszPluscPlusf.boundaryField()[patchI];
                    forAll(DPluszPluscPlusf_, i) { DPluszPluscPlusf_[i] = 0; }
                }
            }

            forAll(cMinus.boundaryField(),patchI)
            {
                if(cMinus.boundaryField()[patchI].type() == "zeroIonicFlux_nonD")//  or "zeroIonicFlux_nonD") // "zeroIonicFlux_nonD" or 
                {
                    //Info << "zeroIonicFlux_nonD ############################"<<endl;
                    scalarField& DMinuszMinusphiE_ = DMinuszMinusphiE.boundaryField()[patchI];
                    forAll(DMinuszMinusphiE_, i) { DMinuszMinusphiE_[i] = 0; }
                    scalarField& DMinuszMinuscMinusf_ = DMinuszMinuscMinusf.boundaryField()[patchI];
                    forAll(DMinuszMinuscMinusf_, i) { DMinuszMinuscMinusf_[i] = 0; }
                }
            }

            if(solveNS)
            {
                if(splitpsi)
                {
                    #include "PNPNSsplit.H"
                }
                else
                {
                    Info << "PNPNS!" << endl;
                    #include "PNPNS.H"
                }
            }
            else
            {
                if(splitpsi)
                {
                    #include "PNPsplit.H"
                }
                else
                {
                    #include "PNP.H"
                }
            }

            forAll(cPlus.boundaryField(),patchI)
            {
                if(cPlus.boundaryField()[patchI].type() == "fixedIonicFlux")
                {
                    scalarField& Cb1 = cPlus.boundaryField()[patchI];
                    const tmp<scalarField>&  Ci1 = cPlus.boundaryField()[patchI].patchInternalField();
                    scalarField& Phib1 = psiE.boundaryField()[patchI];
                    const tmp<scalarField>& Phii1 = psiE.boundaryField()[patchI].patchInternalField();
                    Cb1 = Ci1/(1.0 + zPlus.value()*(Phib1 - Phii1));
                }
            }

            forAll(cMinus.boundaryField(),patchI)
            {
                if(cMinus.boundaryField()[patchI].type() == "fixedIonicFlux")
                {
                    //Info<<nl<<"********************************zero!!!"<<endl;

                    scalarField& Cb2 = cMinus.boundaryField()[patchI];
                    const tmp<scalarField>& Ci2 = cMinus.boundaryField()[patchI].patchInternalField();
                    scalarField& Phib2 = psiE.boundaryField()[patchI];
                    const tmp<scalarField>&  Phii2 = psiE.boundaryField()[patchI].patchInternalField();
                    Cb2 = Ci2/(1.0 + zMinus.value()*(Phib2 - Phii2));
                }
            }
        }
        
        Info<< "maxResidual = " << maxResidual << nl
            << endl;

        // turbulence->correct();
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    } // Time loop closed

    Info<< "End\n" << endl;

    return 0;
}
