//
//AnalysisSort.cpp
//Clustering
//
//Created by Yanina on 07.09.18.
//
#include "AnalysisSort.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include "TBranch.h"
#include "TempArrays.h"
#include <algorithm>

using namespace std;

AnalysisSort::AnalysisSort(TString  InputFilename,TString  OutputFilename) : inTree(0){

}

AnalysisSort::~AnalysisSort(){
}

void AnalysisSort::Loop(TString  InputFilename,TString  OutputFilename) {
    
    // Load Input file and Input Tree //
    TFile *inFile = new TFile(InputFilename, "READ");
    TTree *inTree = (TTree*)inFile->Get("events/events");
    
    // Load Output file and Output Tree //
    TFile *fOutTree = new TFile(OutputFilename,"recreate");
    TTree *outTree = new TTree("events","Clusterized events by Z resolution");
    
    inTree->SetBranchAddress("eventid", &eventid);
    inTree->SetBranchAddress("etot", &etot);
    inTree->SetBranchAddress("nsteps", &nsteps);
    inTree->SetBranchAddress("e_pri", &e_pri);
    inTree->SetBranchAddress("xp", &xp);
    inTree->SetBranchAddress("yp", &yp);
    inTree->SetBranchAddress("zp", &zp);
    inTree->SetBranchAddress("ed", &ed);
    inTree->SetBranchAddress("xp_pri", &xp_pri);
    inTree->SetBranchAddress("yp_pri", &yp_pri);
    inTree->SetBranchAddress("zp_pri", &zp_pri);
    inTree->SetBranchAddress("type_pri", &type_pri);
    
    Long64_t nentries = inTree -> GetEntriesFast();
    
    //output tree variables : set branches
    outTree->Branch("EventId", &EventId, "EventId/I");
    outTree->Branch("n_steps", &n_steps, "n_steps/I");
    outTree->Branch("Xp", &Xp, "Xp/F");
    outTree->Branch("Yp", &Yp, "Yp/F");
    outTree->Branch("Zp", &Zp, "Zp/F");
    outTree->Branch("xpri", &xpri, "xpri/F");
    outTree->Branch("ypri", &ypri, "ypri/F");
    outTree->Branch("zpri", &zpri, "zpri/F");
    outTree->Branch("epri", &epri, "epri/F");
    outTree->Branch("Etot", &Etot, "Etot/F");
    outTree->Branch("typepri", "vector<string>", &typepri);
    
    for (Long64_t jentry = 0; jentry<nentries; jentry++){

        nb = inTree->GetEntry(jentry);

        nbytes += nb;
       
        // Temp Arrays //
        Int_t nTemp = 0;
        const Int_t tempMax = 10000; // if this is not big enough for all the steps there is a risk of overflow.
        
        //PreCluster Arrays //
        vector<TempArrays> PreCluster(tempMax);
        //We have to initialize this vector somehow per event to avoid the problem with the sorting. Just in the coordinate Z 
        for(int p=0;p<tempMax;p++){
        PreCluster[p].Temp_z = 10000.0; //something super big!!!
        }

        // PleCluster Array //
        TempArrays Temp;

        // Final Array //
        vector<TempArrays> FinalArray(tempMax); 
        
        if ( jentry % 1000 ==0 )             cout << jentry << " Entries  " << endl;

        if( nsteps == 0 ) continue;   
        if( etot == 0 ) continue;      
 
        for ( Int_t i=0; i<nsteps; ++i){
            
            if ( isnan((*xp)[i]) ) continue;
            Float_t xn = (*xp)[i];
            Float_t yn = (*yp)[i];
            Float_t zn = (*zp)[i];
            Float_t e_tot = (*ed)[i];
            

            if( e_tot == 0. ) continue; // Only STEPS with some energy deposited //
            if( !InsideFiducialVolume(xn,yn,zn) ) continue; //Only STEPS with energy deposition inside the fiducial volume //
            
            // Fill the Temp arrays with the information in each step //
            PreCluster[nTemp].Temp_x = xn;
            PreCluster[nTemp].Temp_y = yn;
            PreCluster[nTemp].Temp_z = zn;
            PreCluster[nTemp].Temp_etot = e_tot;
            // Move to the next step //
            nTemp++;
            // But if there are lots of depositions, stop counting (just for optimization means)
            if(nTemp >= tempMax) break;   
        }

        if ( nTemp < 1 ) continue;  //avoid events with no good steps

        sort(PreCluster.begin(),PreCluster.end()); // Sorting the temporary array by Z //
        
        // Initialise AfterCluster first value with the PreCluster structure //
        Int_t nOut = 0;
        Temp.Temp_x = PreCluster[nOut].Temp_x*PreCluster[nOut].Temp_etot;
        Temp.Temp_y = PreCluster[nOut].Temp_y*PreCluster[nOut].Temp_etot;
        Temp.Temp_z = PreCluster[nOut].Temp_z*PreCluster[nOut].Temp_etot;
        Temp.Temp_etot = PreCluster[nOut].Temp_etot;
 
        // Initialise z values for steps //       
        z_step = 0;
        zprev_step = 0;
        
        // Now the iteration is only on the events that deposited energy and inside the Fiducial Volume //
        if( nsteps >1){
            for( Int_t step = 1; step<nTemp ; step++){
                z_step = PreCluster[step].Temp_z;
                zprev_step = PreCluster[step-1].Temp_z;
               
                if( MinimumSeparation( z_step, zprev_step)){
                    Temp.Temp_x +=PreCluster[step].Temp_x * PreCluster[step].Temp_etot;
                    Temp.Temp_y +=PreCluster[step].Temp_y * PreCluster[step].Temp_etot;
                    Temp.Temp_z +=PreCluster[step].Temp_z * PreCluster[step].Temp_etot;
                    Temp.Temp_etot +=PreCluster[step].Temp_etot;

                }
                else {
                    // It's a separate cluster, save the previous scatter
                    FinalArray[nOut].Temp_x = Temp.Temp_x/Temp.Temp_etot;
                    FinalArray[nOut].Temp_y = Temp.Temp_y/Temp.Temp_etot;
                    FinalArray[nOut].Temp_z = Temp.Temp_z/Temp.Temp_etot;
                    FinalArray[nOut].Temp_etot = Temp.Temp_etot;
                    nOut++;
                    
                    // And Initialize new cluster
                    Temp.Temp_x =PreCluster[step].Temp_x * PreCluster[step].Temp_etot;
                    Temp.Temp_y =PreCluster[step].Temp_y * PreCluster[step].Temp_etot;
                    Temp.Temp_z =PreCluster[step].Temp_z * PreCluster[step].Temp_etot;
                    Temp.Temp_etot =PreCluster[step].Temp_etot;
                }
            }
        }
        
        // Now we save the last scatter that ocurred //
        if(nTemp>0){
            FinalArray[nOut].Temp_x = Temp.Temp_x/Temp.Temp_etot;
            FinalArray[nOut].Temp_y = Temp.Temp_y/Temp.Temp_etot;
            FinalArray[nOut].Temp_z = Temp.Temp_z/Temp.Temp_etot;
            FinalArray[nOut].Temp_etot = Temp.Temp_etot;
        }
        //Now we finally can fill the new Tree //

        Int_t n_scatters = nOut + 1; 
        Int_t n_steps = 0;
        xpri = xp_pri;
        ypri = yp_pri;
        zpri = zp_pri;
        epri = e_pri;
        EventId = eventid;
        typepri = type_pri;

        // Save in the tree clusterized events //
        //if( n_scatters > 1 ) continue;   //This is to select only single scattering events 
 
        for( Int_t out = 0; out<=nOut ; out++){ 
                Xp = FinalArray[out].Temp_x;
                Yp = FinalArray[out].Temp_y;
                Zp = FinalArray[out].Temp_z;
                Etot = FinalArray[out].Temp_etot ;
                outTree -> Fill();
                n_steps++;
            
        }
    }

    fDebug = 0;
    if (fDebug == 1 )    outTree -> Print();

    outTree -> Write();
    fOutTree -> Close();
}

bool AnalysisSort::MinimumSeparation(Float_t z_pos1, Float_t z_pos2){
    Double_t Separation_Z = 3.; // in mm, first try with 3 mm
    Double_t AverageZ = (z_pos1+z_pos2)/2.0;
    Double_t AvenrageZ_RefToTheBottom = AverageZ + 71.35; // in mm ( O would be the bottom part of the TPC )
    Double_t AvenrageZ_RefToTheBottom_cm = (AverageZ + 71.35)/10.0;  // in cm 
    //Double_t Separation_Z = 3.5 + 0.047*AvenrageZ_RefToTheBottom_cm;  // in mm (for 400 p.e.)
    //Double_t Separation_Z = 4.6 + 0.073*AvenrageZ_RefToTheBottom_cm;  // in mm (for 5000 p.e.)
    
    bool output = false;
    
    if(fabs( z_pos1 - z_pos2) < Separation_Z){
        output = true;
    }
    return output;
}

bool AnalysisSort::InsideFiducialVolume(Float_t x_pos = 0., Float_t y_pos = 0., Float_t z_pos = 0.){
    Float_t SensitiveVolumeRadius = 35.0 ; // In mm, radius sensitive volume, standard tpc size
    Float_t SensitiveVolumeLength = 71.2 ; // In mm, length sensitive volume, standard tpc size
    Float_t Zoffset = -71.35; // In mm, position lower edge sensitve volume
    bool output = false;
    
    if(sqrt(x_pos*x_pos + y_pos* y_pos) < SensitiveVolumeRadius){
          if( z_pos >= Zoffset && z_pos <= SensitiveVolumeLength+Zoffset ) {
            output = true;
        }
    }
    
    return output;
}

Long64_t AnalysisSort::LoadTree(Long64_t entry){
	if(!inTree) return -1;
	Long64_t centry = inTree -> LoadTree(entry);
	if ( inTree -> GetTreeNumber() != fCurrent){
		fCurrent = inTree -> GetTreeNumber();
	}
	return centry;
}
