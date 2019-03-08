//
//  Clustering.cpp
//  Clustering
//
//  Created by Yanina on 07.09.18.
//
//C++ Header Files
#include <string>
#include <sstream>
#include <unistd.h>

#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <TROOT.h>
#include <iostream>
#include <fstream>

#include "AnalysisSort.h"

using namespace std;

TString InputFilename = 'f';
TString Working_Dir = 'f';

void usage();
int main(int argc, char **argv)
{
    int c = 0;
    
    while((c = getopt(argc,argv,"d:i:d:i")) != -1)
    {
        switch(c)    {
            case 'i':
                InputFilename = optarg;
                break;
            case 'd':
                Working_Dir = optarg;
                break;

        }
    }
    
    // Control that the commandline arguments were correctly parsed //
    if(InputFilename=='f' | Working_Dir =='f' ){
        cout << "Command line arguments not correctly parsed!" << endl;
        usage();
    }
    
    //modified by Patricia to run locally
    TString LocalDir = "/home/ab602/Thesis/Simulation_Data/Calibration_simulations_pointsource/";
    TString long_InputFile = Working_Dir + InputFilename + ".root";
    //TString OutputDir = Working_Dir + "proc" + "/"; //modified by Patricia to run locally
    TString OutputDir = LocalDir + "proc" + "/";
    TString OutputRoot = OutputDir + InputFilename + "_Sorted.root";

    cout << "  Your Output Directory will be      " <<  OutputDir << endl;
    
    cout << OutputRoot << "  will be the processed file from  " << long_InputFile <<  "    Dataset" << endl;

    AnalysisSort MyAnalysis(long_InputFile, OutputRoot);
    MyAnalysis.Loop(long_InputFile, OutputRoot);
    
    return 0;

}

void
usage()
{
    exit(0);
}
