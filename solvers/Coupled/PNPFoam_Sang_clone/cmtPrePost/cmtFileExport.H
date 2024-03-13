#ifndef CMTFILEEXPORT_H_INCLUDED
#define CMTFILEEXPORT_H_INCLUDED
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <queue>
#include <stdarg.h>
using namespace std;
void AddDataToFile(std::string fname,std::string header,bool clear_old_file, int nEntries, ...);
void AddDataToFile(std::string fname,std::string header,double rowEntry[],int nEntries,bool clear_old_file=false);


#endif
