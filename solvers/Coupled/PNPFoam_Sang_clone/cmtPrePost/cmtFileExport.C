#include "cmtFileExport.H"

void AddDataToFile(string fname,string header,bool clear_old_file, int nEntries, ...)
{
    bool empty_file=true;
    ifstream ifs(fname.c_str());
    if (ifs.is_open())
    {
      string firstline;
      getline(ifs,firstline);
      if(firstline.compare(header)==0) empty_file=false;
      ifs.close();
    }

    ofstream outFile;
    if(empty_file or clear_old_file)
      outFile.open(fname.c_str(),ios_base::ate);
    else
      outFile.open(fname.c_str(),ios_base::app);

    if(empty_file or clear_old_file)
    {
      outFile<<header<<endl;
    }

    va_list ap;
    va_start(ap, nEntries); //Requires the last fixed parameter (to get the address)
    for(int j=0; j<nEntries; j++)
    {
      if(j<nEntries-1)
        outFile<<setprecision(10)<<scientific<<va_arg(ap,double)<<",";
      else
        outFile<<setprecision(10)<<scientific<<va_arg(ap,double)<<endl;
    }
    va_end(ap);
    outFile.close();
}
void AddDataToFile(std::string fname,std::string header,double rowEntry[],int nEntries,bool clear_old_file)
{
    bool empty_file=true;
    ifstream ifs(fname.c_str());
    if (ifs.is_open())
    {
      string firstline;
      getline(ifs,firstline);
      if(firstline.compare(header)==0) empty_file=false;
      ifs.close();
    }

    ofstream fluxes_file;
    if(empty_file or clear_old_file)
      fluxes_file.open(fname.c_str(),ios_base::ate);
    else
      fluxes_file.open(fname.c_str(),ios_base::app);

    if(empty_file or clear_old_file)
    {
      fluxes_file<<header<<endl;
      for (int j=0;j<nEntries;j++)
      {
          if(j<nEntries-1)
          fluxes_file<<setprecision(10)<<scientific<<rowEntry[j]<<",";
          else fluxes_file<<setprecision(10)<<scientific<<rowEntry[j]<<endl;
      }
    }else
    {
      for (int j=0;j<nEntries;j++)
      {
          if(j<nEntries-1)
          fluxes_file<<setprecision(10)<<scientific<<rowEntry[j]<<",";
          else fluxes_file<<setprecision(10)<<scientific<<rowEntry[j]<<endl;
      }
    }
    fluxes_file.close();
}
