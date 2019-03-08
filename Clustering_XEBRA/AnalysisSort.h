//
//  AnalysisSort.h
//  Clustering
//
//  Created by Yanina on 07.09.18.
//

#ifndef AnalysisSort_h
#define AnalysisSort_h

#include <stdio.h>
#include <vector>
#include "TempArrays.h"
#include <TTree.h>
#include <stddef.h>
#include <string>

using namespace std;

class AnalysisSort{
public:
    // Constructor & Destructor //
    AnalysisSort(TString  InputFilename,TString  OutputFilename);
    ~AnalysisSort();
    
    // Debug //
	Int_t fDebug;
    
    // Control over Tree id //
    virtual Long64_t LoadTree(Long64_t entry);
    Int_t fCurrent;

    Long64_t nbytes = 0, nb = 0;

    // Input Tree //
    static TString InputFilename;
    TTree *inTree;
    vector<float>   *xp = 0;
    vector<float>   *yp = 0;
    vector<float>   *zp = 0;
    vector<float>   *ed = 0;
    Int_t           eventid;
    Float_t         e_pri;
    Int_t           ntpmthits;
    Int_t           nbpmthits;
    Int_t           nLSpmthits;
    Int_t           nWaterpmthits;
    Float_t         etot;
    Int_t           nsteps;
    Float_t         xp_pri;
    Float_t         yp_pri;
    Float_t         zp_pri;
    vector<string>  *type_pri;
    
    // Output Tree //
    TFile *fOutTree;
    static TString OutputFilename;
    TTree *outTree;
    Int_t           EventId;
    Int_t           n_steps;
    Float_t         Xp;
    Float_t         Yp;
    Float_t         Zp;
    Float_t         Etot;
    Float_t         zpri;
    Float_t         xpri;
    Float_t         ypri;
    Float_t         epri;
    vector<string>  *typepri;

    // Loop over events //
    void Loop(TString InputFilename,TString  OutputFilename);
    
    // Cluster by distance in Z //
    bool MinimumSeparation(Float_t z_pos1, Float_t z_pos2);
    Double_t Separation_Z;
    Float_t z_step;
    Float_t zprev_step;
    
    // Temp arrays ( select, sort and cluster phases) //
    Int_t tempMax;
    TempArrays Temp;
    
    // Select fiducial volume //
    bool InsideFiducialVolume(Float_t x_pos, Float_t y_pos, Float_t z_pos);
    Float_t SensitiveVolumeRadius; // In mm
    Float_t SensitiveVolumeLength; // In mm
    Float_t Zoffset; // In mm
};





#endif /* AnalysisSort_h */
