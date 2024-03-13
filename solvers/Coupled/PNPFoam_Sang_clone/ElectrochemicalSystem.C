#include "ElectrochemicalSystem.H"
#include "wordList.H"
#include "IFstream.H"
#include <iostream>
#include <iomanip>
namespace Foam
{
	defineTypeNameAndDebug(ElectrochemicalSystem, 0);
}

template<class GeoField>
void Foam::ElectrochemicalSystem::print
(
	const char* msg,
	Ostream& os,
	const PtrList<GeoField>& flds
)
{
    if (flds.size())
    {
        os << msg;
        forAll(flds, i)
        {
            os<< ' ' << flds[i].name();
        }
        os << endl;
    }
}

void Foam::ElectrochemicalSystem::print(Ostream& os, const wordList& flds)
{
    forAll(flds, i)
    {
        os<< ' ' << flds[i];
    }
    os << endl;
}

Foam::ElectrochemicalSystem::ElectrochemicalSystem
(
	const fvMesh& mesh,
	const dictionary& dict
):
mesh_(mesh),
objects(mesh_, mesh_.time().timeName()),
electricdict_(dict),
phiStart_(readScalar(dict.subDict("phiDyMBoundary").lookup("phiStart"))),
phiEnd_(readScalar(dict.subDict("phiDyMBoundary").lookup("phiEnd"))),
phiInterval_(readScalar(dict.subDict("phiDyMBoundary").lookup("phiInterval"))),
phiInstant_(phiStart_),
writeHeader_(true),
solveNS_(false),
l0(dict.subDict("Electrolyte").lookup("l0")),
C0(dict.subDict("Electrolyte").lookup("C0")),
T(dict.subDict("Electrolyte").lookup("T")),
e(dict.subDict("Electrolyte").lookup("e")),
kB(dict.subDict("Electrolyte").lookup("kB")),
F(dict.subDict("Electrolyte").lookup("F")),
NA(dict.subDict("Electrolyte").lookup("NA")),
eps0(dict.subDict("Electrolyte").lookup("eps0")),
epsr(dict.subDict("Electrolyte").lookup("epsr")),
wsc(dict.subDict("Electrolyte").lookup("wsc")),
rho(dict.subDict("Electrolyte").lookup("rho")),
mu(dict.subDict("Electrolyte").lookup("mu")),
nu(dict.subDict("Electrolyte").lookup("nu")),
z1(dict.subDict("Ions").subDict("C1").lookup("z")),
D1(dict.subDict("Ions").subDict("C1").lookup("D")),
z2(dict.subDict("Ions").subDict("C2").lookup("z")),
D2(dict.subDict("Ions").subDict("C2").lookup("D")),
U_(electricdict_.subDict("Electrolyte").lookup("Velocity")), ///*sangpv add
P_(electricdict_.subDict("Electrolyte").lookup("Pressure")), ///*sangpv add
C2_(electricdict_.subDict("Electrolyte").lookup("Anion")), ///*sangpv add
C1_(electricdict_.subDict("Electrolyte").lookup("Cation")), ///*sangpv add
Phi_(electricdict_.subDict("Electrolyte").lookup("Potential")), ///*sangpv add
Mem1Name(electricdict_.subDict("Desalination").lookup("Mem1Name")), ///*sangpv add
Mem2Name(electricdict_.subDict("Desalination").lookup("Mem2Name")), ///*sangpv add
desaltType(electricdict_.subDict("Desalination").lookup("desaltType")) ///*sangpv add
{

Info<< "ElectricProperties" << endl;
Info << "****************************************************************" << endl;



Info << "length-scale: " << l0.value() << endl;
Info << "Bulk Concentration :" << C0.value() << endl;
Info << "Absolute temperature: " << T.value() << endl;
Info << "elementary charge: " << e.value() << endl;
Info << "Bolt-Zmann constant: " << kB.value() << endl;
Info << "Faraday constant: " << F.value() << endl;
Info << "Avogadro's number: " << NA.value() << endl;
Info << "Vacumm permittivity: " << eps0.value() << endl;
Info << "Dielectric constant: " << epsr.value() << endl;
Info << "Wall Surface Charge: " << wsc.value() << endl;
Info << "Density: " << rho.value() << endl;
Info << "Dynamics viscosity: " << mu.value() << endl;
Info << "Kinematic viscosity: " << nu.value() << endl;
Info << "Diffusivity D1: " << D1.value() << endl;
Info << "Diffusivity D2: " << D2.value() << endl;

    scalar D0 = 0.5*(D1.value() + D2.value());
    DimlD1 = D1.value()/D0;
    DimlD2 = D2.value()/D0;

    scalar Phi0 = (kB.value()*T.value())/(z1.value()*e.value());
    scalar U0 = epsr.value()*eps0.value()*sqr(Phi0)/(mu.value()*l0.value());
    Pe = U0*l0.value()/D0;
    desalinationInletsList =  wordList(dict.subDict("Desalination").lookup("desalinationInlet"));
    desalinationOutletsList = wordList(dict.subDict("Desalination").lookup("desalinationOutlet"));

}
void Foam::ElectrochemicalSystem::saveDesaltPhysics()
{
    scalar v1=0,v2=0;
    scalar flowrate=0,ce=0,rr=0,epir=0,pc=0;
    scalar c_in = 0.0, qc_in = 0.0,c_out = 0.0, qc_out = 0.0,q_in = 0,q_out = 0;

    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (mesh_.boundaryMesh()[patchI].name() == Mem1Name)
        {
            label patchID = patchI;
            const volScalarField& Phi = mesh_.lookupObject<volScalarField>(Phi_);
            const scalarField& PhiBound = Phi.boundaryField()[patchID];
            forAll(PhiBound, i)
            {
                v1 = PhiBound[i];
                break;
            }
        }
        if (mesh_.boundaryMesh()[patchI].name() == Mem2Name)
        {
            label patchID = patchI;
            const volScalarField& Phi = mesh_.lookupObject<volScalarField>(Phi_);
            const scalarField& PhiBound = Phi.boundaryField()[patchID];
            forAll(PhiBound, i)
            {
                v2 = PhiBound[i];
                break;
            }
        }
    }

    scalar I_BottomMem_C1 = I1(Mem1Name);
    scalar I_BottomMem_C2 = I2(Mem1Name);
    scalar I_BottomMem_total = I_BottomMem_C1 + I_BottomMem_C2;
    scalar I_TopMem_C1 = I1(Mem2Name);
    scalar I_TopMem_C2 = I2(Mem2Name);//??
    scalar I_TopMem_total = I_TopMem_C1 + I_TopMem_C2;
    scalar I_Inlet_C1 = I(desalinationInletsList,C1_);
    scalar I_Inlet_C2 = I(desalinationInletsList,C2_);
    scalar I_Inlet_total = I_Inlet_C1 + I_Inlet_C2;
    scalar I_Outlet_C1 = I(desalinationOutletsList,C1_);
    scalar I_Outlet_C2 = I(desalinationOutletsList,C2_);
    scalar I_Outlet_total = I_Outlet_C1 + I_Outlet_C2;
    calcC(C1_, qc_in, q_in, desalinationInletsList);
    calcC(C1_, qc_out, q_out, desalinationOutletsList);

    forAll(desalinationInletsList, patchJ)
    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (mesh_.boundary()[patchI].name() == desalinationInletsList[patchJ])
        {
            scalar nFacesPatchI = mesh_.boundary()[patchI].Cf().size();
            for(int fI = 0; fI < nFacesPatchI; fI++)
            {
                flowrate += massflux(patchI, fI);
            }
        }
    }

    c_in = qc_in/Pe/q_in;
    c_out = qc_out/Pe/q_out;
    if(desaltType == "ED")
    {
      ce = (mag(qc_in) - mag(qc_out))/mag(I_BottomMem_C1);
    }
    else
    {
      ce = 2.0*(mag(qc_in) - mag(qc_out))/mag(I_BottomMem_C1);
    }
    rr =(mag(qc_in) - mag(qc_out))/mag(qc_in);
    epir = (v2 -v1)/ce*2.0;
    pc = (v2 - v1)*I_BottomMem_C1;

    Info<< "Save desaltPhysics..." <<endl;

    fileName outputDir(fileName::null);
    outputDir = mesh_.time().path()/"physicData";
    mkDir (outputDir);
    //-Create file

    const fileName outputFile(outputDir/"desaltPhysics.csv");

    Info << "{! DESALNINATION INFO: flowrate : " << flowrate <<" --- qc_in : " << qc_in << " --- qc_out : " << qc_out << " --- I_BottomMem_C1 : " << I_BottomMem_C1 << endl;
    Info << "                       CE : " << ce <<" --- SRR : " << rr << " --- EPIR : " << epir << " --- PC : " << pc << "}"<<endl;
    if (writeHeader_)
    {
      OFstream os(outputFile,std::ofstream::out | ios_base::ate);
      os << "Voltage,FlowRate,CE,SRR,EPIR,PC,C_in,C_out,I_"<<Mem1Name<<"_C1,I_"<<Mem1Name<<"_C2,I_BottomMem_total,I_"<<Mem2Name<<"_C1,I_"<<Mem2Name<<"_C2,I_TopMem_total,I_Inlet_C1,I_Inlet_C2,I_Inlet_total,I_Outlet_C1,I_Outlet_C2,I_Outlet_total"<<"\r\n";
      os  << v2 - v1 << "," << mag(flowrate) << "," << ce << "," << rr << "," << epir << "," << pc << "," << mag(c_in) << "," << mag(c_out)
        << "," << I_BottomMem_C1 << ","<< I_BottomMem_C2 << ","<< I_BottomMem_total
        << "," << I_TopMem_C1 << "," << I_TopMem_C2 << "," << I_TopMem_total
        << "," << I_Inlet_C1 << "," << I_Inlet_C2 << "," << I_Inlet_total
        << "," << I_Outlet_C1 << "," << I_Outlet_C2 << "," << I_Outlet_total<< "\r\n";
    }
    else
    {
      OFstream os(outputFile,std::ofstream::out | std::ofstream::app);
      os  << v2 - v1 << "," << mag(flowrate) << "," << ce << "," << rr << "," << epir << "," << pc << "," << mag(c_in) << "," << mag(c_out)
        << "," << I_BottomMem_C1 << ","<< I_BottomMem_C2 << ","<< I_BottomMem_total
        << "," << I_TopMem_C1 << "," << I_TopMem_C2 << "," << I_TopMem_total
        << "," << I_Inlet_C1 << "," << I_Inlet_C2 << "," << I_Inlet_total
        << "," << I_Outlet_C1 << "," << I_Outlet_C2 << "," << I_Outlet_total<< "\r\n";
    }
    writeHeader_ = false;
}
void Foam::ElectrochemicalSystem::saveIV()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    word patchName;
    word header = "Voltage";
    forAll(mesh_.boundary(), patchI)
    {
      header += ", " + mesh_.boundary()[patchI].name();
    }
    cout<<" "<<header;
    scalar I_bound[mesh_.boundary().size()];
    forAll(mesh_.boundary(), patchI)
    {
      patchName = mesh_.boundary()[patchI].name();
      I_bound[patchI] = I1(patchName) + I2(patchName);
      Info<<patchName <<" "<<I_bound[patchI]<<endl;
    }
    fileName outputDir(fileName::null);
    outputDir = mesh_.time().path()/"physicData";
    mkDir (outputDir);

    const fileName outFile(outputDir/"IVresponse.csv");

    word vl;
    IFstream is(outFile);
    if(!is.fail())    is>>vl;

    if(vl!="Voltage,")
    {
      OFstream os(outFile,std::ofstream::out);
      os<<header<<"\r\n";
      os<< phiInstant();
      for(int patchI = 0;patchI<mesh_.boundary().size();patchI++)
        os<< ", \t"<< I_bound[patchI];
      os<< "\r\n";
    }else
    {
      OFstream os(outFile,std::ofstream::out | std::ofstream::app);
      os<< phiInstant();
      for(int patchI = 0;patchI<mesh_.boundary().size();patchI++)
        os<< ", \t"<< I_bound[patchI];
      os<< "\r\n";
    }
}


void Foam::ElectrochemicalSystem::changePhi()
{
    phiInstant_ += phiInterval_;
}

void Foam::ElectrochemicalSystem::saveDataVTK(const fvMesh& mesh, int timestep)
{
    Info<< "Save dataVTK..." <<endl;

    //- Create directory
    fileName outputDir(fileName::null);
    outputDir = mesh.time().path()/"dataVTK";
    mkDir (outputDir);

    // Search for list of objects for this time
    //IOobjectList objects(mesh_, mesh_.time().timeName());
    // Construct the vol fields and point fields
    HashSet<word> selectedFields;
    vtkMesh vMesh(const_cast<fvMesh&>(mesh));

    PtrList<volScalarField> vsf;
    PtrList<volVectorField> vvf;
    PtrList<pointScalarField> psf;
    PtrList<pointVectorField> pvf;
    readFields(vMesh, mesh, objects, selectedFields, vsf);
    readFields(vMesh, mesh, objects, selectedFields, vvf);
    readFields(vMesh, pointMesh::New(mesh), objects, selectedFields, psf);
    readFields(vMesh, pointMesh::New(mesh), objects, selectedFields, pvf);


    print("    volScalarFields            :", Info, vsf);
    print("    volVectorFields            :", Info, vvf);
    print("    pointScalarFields          :", Info, psf);
    print("    pointVectorFields          :", Info, pvf);

    label nVolFields   = vsf.size() + vvf.size();
    label nPointFields = psf.size() + pvf.size();

    //- Create file vtk
    string vtkName = mesh.time().caseName();
    std::ostringstream os,str1; os << phiInstant_ ;
    string phiName = os.str();
    word timename = Foam::name(timestep-1);

    fileName vtkFileName;
    if( timestep<2)
      vtkFileName = (outputDir/vtkName + "_Phi_" + phiName + ".vtk");
    else vtkFileName = (outputDir/vtkName + "_Phi_" + phiName + "_step_"+timename+".vtk");

    //- write mesh
    internalWriter writer(vMesh, false, vtkFileName);

    //- write cell data
    writeFuns::writeCellDataHeader
    (
        writer.os(),
        vMesh.nFieldCells(),
        1+nVolFields
    );
    writer.writeCellIDs();
    writer.write(vsf);
    writer.write(vvf);

    //- write point data
    writeFuns::writePointDataHeader
    (
        writer.os(),
        vMesh.nFieldPoints(),
        nVolFields+nPointFields
    );

    // pointFields
    writer.write(psf);
    writer.write(pvf);

    // Interpolated volFields
    volPointInterpolation pInterp(mesh_);
    writer.write(pInterp, vsf);
    writer.write(pInterp, vvf);
}

void Foam::ElectrochemicalSystem::saveDataVTK(string filnename, PtrList<volScalarField> vsf)
{
    Info<< "Save dataVTK..." <<endl;
    //- Create directory
    fileName outputDir(fileName::null);
    outputDir = mesh_.time().path()/"dataVTK";
    mkDir (outputDir);
    print("    volScalarFields            :", Info, vsf);

    label nVolFields   = vsf.size();
    ///- Create file vtk
    string vtkName = mesh_.time().caseName();
    std::ostringstream os; os << phiInstant_ ;
    string phiName = os.str();

    fileName vtkFileName(outputDir/vtkName + filnename + ".vtk");
    Info<< "    File name  : " << vtkFileName << endl;

    vtkMesh vMesh(const_cast<fvMesh&>(mesh_));
    //- write mesh
    internalWriter writer(vMesh, false, vtkFileName);

    //- write cell data
    writeFuns::writeCellDataHeader
    (
        writer.os(),
        vMesh.nFieldCells(),
        1+nVolFields
    );
    writer.writeCellIDs();
    writer.write(vsf);
    // Interpolated volFields
    volPointInterpolation pInterp(mesh_);
    writer.write(pInterp, vsf);
}

void Foam::ElectrochemicalSystem::writeIonicFluxBoundary
(
    word fieldName,
    scalar D,
    scalar z,
    scalar phiInstant,
    word patchName
)
{
    Info<< "Save writeIonicFluxBoundary..." <<endl;

    fileName outputDir(fileName::null);

    outputDir = mesh_.time().path()/"physicData";

    mkDir (outputDir);



    label patchID = mesh_.boundaryMesh().findPatchID(patchName);
    vectorField coordinates(mesh_.Cf().boundaryField()[patchID].size());
    vectorField n(mesh_.Cf().boundaryField()[patchID].size());
    forAll(mesh_.Cf().boundaryField()[patchID],i)
    {
        coordinates[i] = mesh_.Cf().boundaryField()[patchID][i];
        n[i]= mesh_.Sf().boundaryField()[patchID][i]/mesh_.magSf().boundaryField()[patchID][i];
    }
    // // add J
    // word fieldName("Phi");
    const volScalarField& field = mesh_.lookupObject<volScalarField>(fieldName);
    const volScalarField& Phi = mesh_.lookupObject<volScalarField>("Phi");
    surfaceScalarField magSf = mag(mesh_.Sf());
    const scalarField& fieldBound = field.boundaryField()[patchID];
    const scalarField fieldInner = field.boundaryField()[patchID].patchInternalField();
    const scalarField& PhiBound = Phi.boundaryField()[patchID];
    const scalarField PhiInner = Phi.boundaryField()[patchID].patchInternalField();
    scalarField& magS = magSf.boundaryField()[patchID];

    const labelList& cellPatch = mesh_.boundaryMesh()[patchID].faceCells();
    vectorField J_patch(mesh_.Cf().boundaryField()[patchID].size());
    forAll(cellPatch, i)
    {
        vector delta = mesh_.C()[cellPatch[i]] - mesh_.Cf().boundaryField()[patchID][i];
        scalar gradC2 = (fieldBound[i] - fieldInner[i])*magS[i]/mag(delta);
        scalar gradPhi = (PhiBound[i] - PhiInner[i])*magS[i]/mag(delta);
        J_patch[i] = -(D*(gradC2 + z*fieldBound[i]*gradPhi)*magS[i] + Pe*fieldBound[i]*massflux(patchID, i))*n[i];
    }

    const fileName outputFile(outputDir/"ionicFlux-" + fieldName + patchName + name(phiInstant) + ".csv");
    OFstream os(outputFile, ios_base::app);
    os  << "faceID, x , y, z, j1x, j1y, j1z" << endl;

    forAll(cellPatch, fI)
    {
        os  << fI << ", "
            << coordinates[fI].x() << ", "
            << coordinates[fI].y() << ", "
            << coordinates[fI].z() << ", "
            << J_patch[fI].x() << ", "
            << J_patch[fI].y() << ", "
            << J_patch[fI].z() << ", "
            <<endl;
    }
}


void Foam::ElectrochemicalSystem::CalPowerDissipation(string filenamehead)
{
  const volScalarField& Phi = mesh_.lookupObject<volScalarField>(Phi_);
  const volScalarField& C1 = mesh_.lookupObject<volScalarField>(C1_);
  const volScalarField& C2 = mesh_.lookupObject<volScalarField>(C2_);
  const volVectorField& U = mesh_.lookupObject<volVectorField>(U_);
  const volScalarField& P = mesh_.lookupObject<volScalarField>(P_);

  volScalarField oD("Ohmic_dissipation",Phi),cD("Spacecharge_dissipation",Phi),vD("Viscous_dissipation",Phi);

  volVectorField gradC1 = fvc::grad(C1);
  volVectorField gradC2 = fvc::grad(C2);
  volVectorField  gradPhi = fvc::grad(Phi);
  volVectorField  gradP = fvc::grad(P);
  volScalarField magGradPhi = mag(gradPhi);
  scalar R(1);
  ///  ohmic dissipation
  forAll(oD.internalField(),i)
  {
    oD.internalField()[i] = pow(F.value(),2)/(R*T.value())*pow(magGradPhi[i],2)*(D1.value()*pow(z1.value(),2)*C1.internalField()[i]+D2.value()*pow(z2.value(),2)*C2.internalField()[i]);
  }

  forAll(cD.internalField(),i)
  {
      ///const label cellCemId = ListCellsOfCEM[icells];
    cD.internalField()[i] =  F.value()*dot(D1.value()*z1.value()*gradC1[i]+D2.value()*z2.value()*gradC2[i],gradPhi[i]);
  }
  ///  viscous dissipation
  forAll(vD.internalField(),i)
  {
    vD.internalField()[i]= U[i]&(gradP[i]+gradPhi[i]*F.value()*(z1.value()*C1.internalField()[i]+z2.value()*C2.internalField()[i]));
  }
  string filnename = "powerDissipation_"+filenamehead;

  PtrList<volScalarField> vsf;

  vtkMesh vMesh(const_cast<fvMesh&>(mesh_));

  addField(vMesh,vMesh,oD,vsf);
  addField(vMesh,vMesh,cD,vsf);
  addField(vMesh,vMesh,vD,vsf);

  saveDataVTK(filnename,vsf);
}
void Foam::ElectrochemicalSystem::CurrentPerTimestep(scalar timestep)
{
    scalar v1=0,v2=0;
    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (mesh_.boundaryMesh()[patchI].name() == Mem1Name)
        {
            label patchID = patchI;
            const volScalarField& Phi = mesh_.lookupObject<volScalarField>(Phi_);
            const scalarField& PhiBound = Phi.boundaryField()[patchID];
            forAll(PhiBound, i)
            {
                v1 = PhiBound[i];
            }
        }
        if (mesh_.boundaryMesh()[patchI].name() == Mem2Name)
        {
            label patchID = patchI;
            const volScalarField& Phi = mesh_.lookupObject<volScalarField>(Phi_);
            const scalarField& PhiBound = Phi.boundaryField()[patchID];
            forAll(PhiBound, i)
            {
                v2 = PhiBound[i];
            }
        }
    }
    scalar I_BottomMem_C1 = I1(Mem1Name);
    scalar I_BottomMem_C2 = I2(Mem1Name);
    scalar I_BottomMem_total = I_BottomMem_C1 + I_BottomMem_C2;

    std::ostringstream os; os << phiInstant_ ;
    string phiName = os.str();
    fileName outputDir(fileName::null);
    outputDir = mesh_.time().path()/"CurrentPerTimestep";
    mkDir (outputDir);

    const fileName outputFile(outputDir/phiName + "V.csv");

    if (timestep==1)
    {
      OFstream os1(outputFile, ios_base::trunc);
      os1 << "dV,timeStep,I,iSum" << endl;
      os1 << v2-v1 << "," << timestep << "," << I_BottomMem_total << "," <<  CurrentSum() << endl;
    }else
    {
      OFstream os1(outputFile, ios_base::app);
      os1 << v2-v1 << "," << timestep<< "," << I_BottomMem_total << "," <<  CurrentSum() << endl;
    }


}
scalar Foam::ElectrochemicalSystem::CurrentPerIteration(scalar i,bool writeHeader)
{
    scalar v1=0,v2=0;
    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (mesh_.boundaryMesh()[patchI].name() == Mem1Name)
        {
            label patchID = patchI;
            const volScalarField& Phi = mesh_.lookupObject<volScalarField>(Phi_);
            const scalarField& PhiBound = Phi.boundaryField()[patchID];
            forAll(PhiBound, i)
            {
                v1 = PhiBound[i];
            }
        }
        if (mesh_.boundaryMesh()[patchI].name() == Mem2Name)
        {
            label patchID = patchI;
            const volScalarField& Phi = mesh_.lookupObject<volScalarField>(Phi_);
            const scalarField& PhiBound = Phi.boundaryField()[patchID];
            forAll(PhiBound, i)
            {
                v2 = PhiBound[i];
            }
        }
    }
    scalar I_BottomMem_C1 = I1(Mem1Name);
    scalar I_BottomMem_C2 = I2(Mem1Name);
    scalar I_BottomMem_total = I_BottomMem_C1 + I_BottomMem_C2;

    std::ostringstream os; os << phiInstant_ ;
    string phiName = os.str();
    fileName outputDir(fileName::null);
    outputDir = mesh_.time().path()/"CurrentPerIteration";
    mkDir (outputDir);

    const fileName outputFile(outputDir/phiName + "V.csv");
    OFstream os1(outputFile, ios_base::app);

    if (writeHeader)
    {
        os1 << "dV,Iteration,I,iSum" << endl;
    }

    os1 << v2-v1 << "," << i+1 << "," << I_BottomMem_total << "," <<  CurrentSum() << endl;
    return I_BottomMem_total;
}

Foam::scalar Foam::ElectrochemicalSystem::CurrentSum()
{
    scalar iSum(0);
    iSum = I1("bottom") + I2("bottom");
    return iSum;
}

Foam::scalar Foam::ElectrochemicalSystem::I(wordList patchNames,word C12)
{
  const polyBoundaryMesh& patches = mesh_.boundaryMesh();
  scalar J12 = 0;
  dimensionedScalar DimlD12(electricdict_.subDict("Ions").subDict(C12).lookup("D"));
  dimensionedScalar z12(electricdict_.subDict("Ions").subDict(C12).lookup("z"));
  forAll(patchNames,nI)
  {
   	J12 += calcIonicFluxJ(patchNames[nI],C12,DimlD12.value(),z12.value());
  }

	return z12.value()*J12;
}

Foam::scalar Foam::ElectrochemicalSystem::I1(word patchName)
{
  bool foundPatch = false;
	forAll(mesh_.boundaryMesh(), patchI)
      if(mesh_.boundary()[patchI].name() == patchName)
      {
        foundPatch = true;
      }
  if(foundPatch)
  {
    scalar J1 = calcIonicFluxJ(patchName,C1_,DimlD1,z1);
    return z1.value()*J1;
  }
  else
  {
    Info<<"No boundary named "<<patchName<<"!"<<endl;
    return 0;
  }
}
Foam::scalar Foam::ElectrochemicalSystem::I2(word patchName)
{
  bool foundPatch = false;
	forAll(mesh_.boundaryMesh(), patchI)
      if(mesh_.boundary()[patchI].name() == patchName)
        foundPatch = true;
  if(foundPatch)
  {
    scalar J2 = calcIonicFluxJ(patchName,C2_,DimlD2,z2);
    return z2.value()*J2;
  }
  else
  {
    Info<<"No boundary named "<<patchName<<"!"<<endl;
    return 0;
  }
}

Foam::scalar Foam::ElectrochemicalSystem::calcIonicFluxJ(word patchName, word fieldName, scalar DimlD, dimensionedScalar z)
{
    scalar diff(0.0), conv(0.0), migr(0.0), sumFlux(0.0);

    const volScalarField& field = mesh_.lookupObject<volScalarField>(fieldName);
	  const volScalarField& Phi = mesh_.lookupObject<volScalarField>(Phi_);
    surfaceScalarField magSf = mag(mesh_.Sf());

	forAll(mesh_.boundaryMesh(), patchI)
	{
    if((mesh_.boundary()[patchI].name() == patchName) && ( mesh_.boundary()[patchI].type() != "empty"))
    {
        label patchID = patchI;
        const scalarField& fieldBound = field.boundaryField()[patchID];
        const scalarField fieldInner = field.boundaryField()[patchID].patchInternalField();
        const scalarField& PhiBound = Phi.boundaryField()[patchID];
        const scalarField PhiInner = Phi.boundaryField()[patchID].patchInternalField();
        scalarField& magS = magSf.boundaryField()[patchID];

        const labelList& cellPatch = mesh_.boundaryMesh()[patchID].faceCells();


      forAll(cellPatch, i)
      {
          vector delta = mesh_.C()[cellPatch[i]] - mesh_.Cf().boundaryField()[patchID][i];
          scalar gradC = (fieldBound[i] - fieldInner[i])/mag(delta);
          scalar gradPhi = (PhiBound[i] - PhiInner[i])/mag(delta);
          diff += -DimlD*gradC*magS[i];
          migr += -DimlD*z.value()*fieldBound[i]*gradPhi*magS[i];
          conv += Pe*fieldBound[i]*massflux(patchID, i);
          scalar flux = -DimlD*(gradC + z.value()*fieldBound[i]*gradPhi)*magS[i] + Pe*fieldBound[i]*massflux(patchID,i);
          sumFlux += flux;
      }
    }
    }

    return sumFlux;
}

void Foam::ElectrochemicalSystem::calcC(word fieldName, scalar& q_con, scalar& q_, wordList physicalpatchName)
{
    const volScalarField& field = mesh_.lookupObject<volScalarField>(fieldName);
    q_con = 0;
    q_ = 0;
    forAll(physicalpatchName, patchJ)
    {
      bool foundBound = false;
      forAll(mesh_.boundaryMesh(), patchI)
      {
          if (mesh_.boundary()[patchI].name() == physicalpatchName[patchJ])
          {
              label patchID = patchI;
              const scalarField& fieldBound = field.boundaryField()[patchID];
              forAll(fieldBound, i)
              {
                  q_con += fieldBound[i]*massflux(patchID, i)*Pe;
                  q_ += massflux(patchID, i);
              }
              foundBound = true;
          }
      }
      if (!foundBound) Info << " ~~~ PNPNSFoam - Warning: No patch found for the patch name "<<physicalpatchName[patchJ]<<". Check the wordList that defined in TransportDict/Desalination!!!"<< endl;
    }
}
Foam::scalar Foam::ElectrochemicalSystem::massflux(label patchID, label faceID)
{
        const surfaceScalarField& phi = mesh_.lookupObject<surfaceScalarField>("phi");
        const scalarField& phiPatch = phi.boundaryField()[patchID];
        return phiPatch[faceID];
}
