
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include "TBranch.h"
#include <algorithm>
#include <math.h>

#include "AnalysisSort.h"
#include "TempArrays.h"

using namespace std;

AnalysisSort::AnalysisSort(){;}  
AnalysisSort::~AnalysisSort(){;} 

void AnalysisSort::Run( vector<float> *xp, vector<float> *yp, vector<float> *zp, vector<float> *ed ) {
       
        // Temp Arrays //
        Int_t nTemp = 0;
        const Int_t tempMax = 500000; // if this is not big enough for all the steps there is a risk of overflow.
        
        //PreCluster Arrays //
        vector<TempArrays> PreCluster(tempMax);
        //We have to inicialize this vector somehow (per event) with a super big number to avoid the problem with the sorting. Just in the coordinate Z 
        for(int p=0;p<tempMax;p++){
        PreCluster[p].Temp_z = 5000.0;
        }

        // PleCluster Array //
        TempArrays Temp;

        // Final Array //
        vector<TempArrays> FinalArray(tempMax); 
       
        for ( Int_t i=0; i<xp->size(); ++i){
            
            if ( isnan((*xp)[i]) ) continue;
            Float_t xn = (*xp)[i];
            Float_t yn = (*yp)[i];
            Float_t zn = (*zp)[i];
            Float_t e_tot = (*ed)[i];
            
            if( e_tot == 0. ) continue;

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

        sort(PreCluster.begin(),PreCluster.end()); // Sorting the temporary array by Z  (minimum to maximum) //

        // Initialise AfterCluster first value with the PreCluster structure //
        Int_t nOut = 0;
        Temp.Temp_x = PreCluster[nOut].Temp_x*PreCluster[nOut].Temp_etot;
        Temp.Temp_y = PreCluster[nOut].Temp_y*PreCluster[nOut].Temp_etot;
        Temp.Temp_z = PreCluster[nOut].Temp_z*PreCluster[nOut].Temp_etot;
        Temp.Temp_etot = PreCluster[nOut].Temp_etot;
 
        // Initialise z values for steps //       
        z_step = 0;
        zprev_step = 0;
        x_step = 0;
        xprev_step = 0;
        y_step = 0;
        yprev_step = 0;
        
        // Now the iteration is only on the events that deposited energy, by default in the LXe Sensitive Volume //
        if( (xp->size()) > 1 ){  //Is this necessary?

            for( Int_t step = 1; step<nTemp ; step++){
                z_step = PreCluster[step].Temp_z;
                zprev_step = PreCluster[step-1].Temp_z;
                x_step = PreCluster[step].Temp_x;
                xprev_step = PreCluster[step-1].Temp_x;
                y_step = PreCluster[step].Temp_y;
                yprev_step = PreCluster[step-1].Temp_y;
               
                if( MinimumSeparation( z_step, zprev_step) and MinimumSeparation_xy(x_step, xprev_step, y_step, yprev_step)){
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

        //At this point we shuold have all the information about the clusters. We should file the vectors Etot_vector

        int n_scatters = nOut + 1; 

       x_vector.clear();
       y_vector.clear();
       z_vector.clear();
       etot_vector.clear(); //clear the vectors before filling just in case

       for( Int_t out = 0; out< n_scatters ; out++){ 
                 x_vector.push_back( FinalArray[out].Temp_x );
                 y_vector.push_back( FinalArray[out].Temp_y );
                 z_vector.push_back(  FinalArray[out].Temp_z);
                 etot_vector.push_back( FinalArray[out].Temp_etot );
        }

} //end of the RUN


vector<float> AnalysisSort::GetXVector(){
  return this->x_vector;
}

vector<float> AnalysisSort::GetYVector(){
  return this->y_vector;
}

vector<float> AnalysisSort::GetZVector(){
  return this->z_vector;
}

vector<float> AnalysisSort::GetEnergyVector(){
  return this->etot_vector;
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

bool AnalysisSort::MinimumSeparation_xy(Float_t x_pos1, Float_t x_pos2, Float_t y_pos1, Float_t y_pos2){
    Double_t Separation_xy = 10.; // in mm 
    //ARBITRARILY CHOOSEN VALUE FOR xy POSITION RECONSTRUCTION ### MUST BE REPLACED BY A PROPER MODEL !!! 
  
    bool output = false;
    
    if( (pow((x_pos1 - x_pos2),2) + pow((y_pos1 - y_pos2),2)) < pow(Separation_xy,2)){
        output = true;
    }
    return output;

}

