#ifndef CMTVTU_H
#define CMTVTU_H

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "fvCFD.H"
#include "wordReList.H"
#include "IOobjectList.H"
#include "readFields.H"

#include "vtkMesh.H"
#include "internalWriter.H"

#define number_width  16
#define number_precision 6
void Save_Data_PVTU_VTUs(const fvMesh& mesh_, std::string filepath,std::string filename,int step=-1);
void Save_Data_VTU(const vtkMesh& vMesh, std::string filepath,int step=-1);
std::string MakeFileName(std::string filename, bool pvtu, int rank, int step=-1);
void WriteDatatypeBegin(ofstream&file,std::string datatype);
void WriteDatatypeEnd(ofstream&file,std::string datatype);
void WriteVTUheader(ofstream &file,const int npoints,const  int ncells, const word gridType);
void WriteVTUFooter(ofstream &file, std::string gridType);
void WritePVTUheader(ofstream &file, std::string gridType);
void WritePVTUFooter(ofstream &file, std::string gridType);
#endif
