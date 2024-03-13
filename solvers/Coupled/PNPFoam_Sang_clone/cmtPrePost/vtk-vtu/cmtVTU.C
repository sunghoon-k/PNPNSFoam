#include "cmtVTU.H"
#include "OFstream.H"


void WriteDatatypeBegin(ofstream&file,string datatype)
{
  file<<"   <"<<datatype<<" format=\"ascii\">"<<endl;
}
void WriteDatatypeEnd(ofstream&file,string datatype)
{
  file<<"   </"<<datatype<<">"<<endl;
}

void WriteVTUheader(ofstream &file,const int npoints,const  int ncells, word gridType)
{
	file<<"<?xml version=\"1.0\"?>"<<endl;
	file<<"<VTKFile type=\"UnstructuredGrid\">"<<endl;
  file<<" <"<<gridType<<">"<<endl;
	file<<"  <Piece NumberOfPoints=\"" <<npoints<<"\" NumberOfCells=\""<<ncells<<"\">"<<endl;
}
void WriteVTUFooter(ofstream &file, string gridType)
{
  file<<"  </Piece>"<<endl;
  file<<" </"<<gridType<<">"<<endl;
  file<<"</VTKFile>"<<endl;
}
void WritePVTUheader(ofstream &file, word gridType)
{
  file<<"<?xml version=\"1.0\"?>"<<endl;
  file<<"<VTKFile type=\"P"<<gridType<<"\" byte_order=\"LittleEndian\">"<<endl;
  file<<" <P"<<gridType<<" GhostLevel=\"0\">"<<endl;
  /// nodes:
  file<<"  <PPoints>"<<endl;
  file<<"   <PDataArray NumberOfComponents=\"3\" type=\"Float32\"  name=\"Position\" format=\"ascii\" />"<<endl;
  file<<"  </PPoints>"<<endl;
  /// cells:
  file<<"  <PCells>"<<endl;
  file<<"   <PDataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\" />"<<endl;
  file<<"   <PDataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\" />"<<endl;
  file<<"   <PDataArray Name=\"types\" type=\"UInt8\" format=\"ascii\" />"<<endl;
  file<<"  </PCells>"<<endl;
}
void WritePVTUFooter(ofstream &file, string gridType)
{
  file<<" </P"<<gridType<<">"<<endl;
  file<<"</VTKFile>"<<endl;
}
string MakeFileName(string filename, bool pvtu, int rank, int step)
{
  string fileName;
  if(pvtu)
  {
    if(step>-1)
      fileName = "MAIN_"+filename + "_step_" + Foam::name(step) + ".pvtu";
    else
      fileName = "MAIN_" + filename + ".pvtu";
  }
  else
  {
   if(step>-1)
      fileName = "Proc_" + Foam::name(rank) + "_"+ filename + "_step_" + Foam::name(step) + ".vtu";
   else
      fileName = "Proc_" + Foam::name(rank) + "_" +filename + ".vtu";
  }
  return fileName;
}
void Save_Data_VTU(const vtkMesh& vMesh,string filePathName,int step)
{
  const fvMesh& mesh = vMesh.mesh();
  const vtkTopo& topo = vMesh.topo();
  const labelList& addPointCellLabels = topo.addPointCellLabels();
  const label nTotPoints = mesh.nPoints() + addPointCellLabels.size();
  const labelListList& vtkVertLabels = topo.vertLabels();

  // Search for list of objects for this time
  IOobjectList objects(vMesh, vMesh.mesh().time().timeName());
  HashSet<word> selectedFields;
  //- write mesh
  internalWriter writer(vMesh, false, filePathName);
    word us("UnstructuredGrid");
  WriteVTUheader(writer.os(),nTotPoints,vMesh.nFieldCells(),us);

  ///WritePointsCoodinate(file,M->Node);:

  DynamicList<floatScalar> ptField(3*nTotPoints);
  writeFuns::insert(mesh.points(), ptField);
  const pointField& ctrs = mesh.cellCentres();
  forAll(addPointCellLabels, api)
  {
      writeFuns::insert(ctrs[addPointCellLabels[api]], ptField);
  }
  writeFuns::write(writer.os(), false , ptField);

  ///WriteCellGeometry(file,M->nCells,M->iNodesofCell,M->NodesofCell,M->GeometricCellType,NULL);

  const labelList& vtkCellTypes = topo.cellTypes();

  writer.os() << "CELL_TYPES " << vtkCellTypes.size() << std::endl;
  // Make copy since writing might swap stuff.
  DynamicList<label> cellTypes(vtkCellTypes.size());
  writeFuns::insert(vtkCellTypes, cellTypes);
  writeFuns::write(writer.os(), false, cellTypes);



  PtrList<volScalarField> vsf;
  PtrList<volVectorField> vvf;
  PtrList<pointScalarField> psf;
  PtrList<pointVectorField> pvf;

  readFields(vMesh, vMesh, objects, selectedFields, vsf);
  readFields(vMesh, vMesh, objects, selectedFields, vvf);
  readFields(vMesh, pointMesh::New(vMesh), objects, selectedFields, psf);
  readFields(vMesh, pointMesh::New(vMesh), objects, selectedFields, pvf);
  ///WriteDataScalars(file,scalarsname,npoints,nscalars,scalars);  if(vec!=NULL)  WriteDataVector(file,vecname,*vec);

  //- write cell data
  string celldata("CellData");
  string pointdata("PointData");
  WriteDatatypeBegin(writer.os(),celldata);
  if (!vsf.empty())
    writer.write(vsf,true);
  if (!vsf.empty())
    writer.write(vsf,true);
  WriteDatatypeEnd(writer.os(),celldata);
  // pointFields
  WriteDatatypeBegin(writer.os(),pointdata);
  if (!psf.empty())
    writer.write(psf,true);
  if (!pvf.empty())
    writer.write(pvf,true);
  // Interpolated volFields
  volPointInterpolation pInterp(vMesh.mesh());
  if (!vsf.empty())
    writer.write(pInterp, vsf,true);
  if (!vvf.empty())
    writer.write(pInterp, vvf,true);
  WriteDatatypeEnd(writer.os(),pointdata);

  WriteVTUFooter(writer.os(),us);
  return;
}
void Save_Data_PVTU_VTUs(const fvMesh& mesh_, string filepath,string filename,int step)
{
  string filePathName,fileName;
  int Rank=Pstream::master();
  int nProcesses;
  vtkMesh vMesh(const_cast<fvMesh&>(mesh_));
  mkDir(filePathName);
  if (Rank==0 and nProcesses>1)
  {
    fileName = MakeFileName(filename,true,Rank,step);
    filePathName = filepath +"/" + fileName;
    ofstream file(filePathName.c_str());
    word us("UnstructuredGrid");
    WritePVTUheader(file,us);
///---
    file<<"<PPointData>"<<endl;

    // Search for list of objects for this time
    IOobjectList objects(mesh_, mesh_.time().timeName());

    // Construct the vol fields and point fields
    HashSet<word> selectedFields;

    PtrList<volScalarField> vsf;
    PtrList<volVectorField> vvf;
    PtrList<pointScalarField> psf;
    PtrList<pointVectorField> pvf;

    readFields(vMesh, vMesh, objects, selectedFields, vsf);
    readFields(vMesh, vMesh, objects, selectedFields, vvf);
    readFields(vMesh, pointMesh::New(vMesh), objects, selectedFields, psf);
    readFields(vMesh, pointMesh::New(vMesh), objects, selectedFields, pvf);

    if (!vsf.empty())
      forAll(vsf, fieldi)    file<<"<PDataArray Name=\""<<vsf[fieldi].name()<<"\" type=\"Float32\" format=\"ascii\" NumberOfComponents=\"1\"/>"<<endl;
    if (!vvf.empty())
      forAll(vvf, fieldi)    file<<"<PDataArray Name=\""<<vvf[fieldi].name()<<"\" type=\"Float32\" format=\"ascii\" NumberOfComponents=\"3\"/>"<<endl;
    if (!psf.empty())
      forAll(psf, fieldi)    file<<"<PDataArray Name=\""<<psf[fieldi].name()<<"\" type=\"Float32\" format=\"ascii\" NumberOfComponents=\"1\"/>"<<endl;
    if (!pvf.empty())
      forAll(pvf, fieldi)    file<<"<PDataArray Name=\""<<pvf[fieldi].name()<<"\" type=\"Float32\" format=\"ascii\" NumberOfComponents=\"3\"/>"<<endl;

    file<<"</PPointData>"<<endl;
///---
    if(step>-1)
    {
    for (int p=0;p<nProcesses;++p)
      file<<"<Piece Source=\""<<filename<<"_proc_"<<p<<"_step_"<<step<<".vtu\" />"<<endl;
    }
    else
    {
    for (int p=0;p<nProcesses;++p)
      file<<"<Piece Source=\""<<filename<<"_proc_"<<p<<".vtu\" />"<<endl;
    }
    WritePVTUFooter(file,us);
    file.close();
  } else
  fileName = MakeFileName(filename,false,Rank,step);

  filePathName = filepath +"/" + fileName;
  /// Save file on all processes
  Save_Data_VTU(vMesh,filePathName,step);
	return;
}
