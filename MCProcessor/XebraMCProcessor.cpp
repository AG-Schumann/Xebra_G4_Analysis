//C++ Header Files
#include <string>
#include <sstream>
#include <unistd.h>
//Root Header Files
#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <TROOT.h>
#include <iostream>
#include <fstream>
//Processing Header
#include "MainProcessor.h"
#include "AnalysisSort.h"

using namespace std;

TString InputFilename = 'f';
TString Working_Dir = 'f';
TString USE_NEST ='f';

void usage();
int main(int argc, char **argv)
{
    int c = 0;
    
    while((c = getopt(argc,argv,"d:i:d:i:n")) != -1)
    {
        switch(c)    {
            case 'i':
                InputFilename = optarg;
                break;
            case 'd':
                Working_Dir = optarg;
                break;
            case 'n':
                USE_NEST = optarg;
                break;
        }
    }
    
    // Control that the commandline arguments were correctly parsed //
    if(InputFilename=='f' | Working_Dir =='f' ){
        cout << "Command line arguments not correctly parsed!" << endl;
        usage();
    }
    TString LocalDir = "/home/alex/Thesis-Copy/Xebra_G4_Analysis/MCProcessor/"; 
    TString long_InputFile = Working_Dir + InputFilename + ".root";
    TString OutputDir = LocalDir + "proc" + "/";
    TString OutputRoot = OutputDir + InputFilename + "_Proc.root";
 
    cout << "  Your Output Directory will be      " <<  OutputDir << endl;
    
    cout << OutputRoot << "  will be the processed file from  " << long_InputFile <<  "    Dataset" << endl;

    // Instance of the main processor
    MainProcessor MyAnalysis(long_InputFile, OutputRoot);
    // Loop over and cluster events
    MyAnalysis.Loop(long_InputFile, OutputRoot);

    
    return 0;

}

void
usage()
{
    exit(0);
}
