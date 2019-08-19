//
// AnalysisSort.h
// Clustering
//
// Created by Yanina on 07.09.18.
//

#ifndef AnalysisSort_h
#define AnalysisSort_h

#include <stdio.h>
#include <vector>
#include <TTree.h>
#include <stddef.h>
#include <string>
#include "TempArrays.h"
#include <math.h>

using namespace std;

class AnalysisSort{
public:
    // Constructor & Destructor //
    AnalysisSort();
   ~AnalysisSort();

    //Functions wich provide a vector as an output //
    vector<float> GetEnergyVector();
    vector<float> GetXVector();
    vector<float> GetYVector();
    vector<float> GetZVector();

    //Intermediate vectors //
    vector<float> x_vector;
    vector<float> y_vector;
    vector<float> z_vector;
    vector<float> etot_vector;

    //variables//
    float z_step;
    float zprev_step;
    
    // Loop over events //
    void Run( vector<float> *xp, vector<float> *yp, vector<float> *zp, vector<float> *ed );
    
    // Cluster by distance in Z //
    bool MinimumSeparation(Float_t z_pos1, Float_t z_pos2);
    Double_t Separation_Z;

    // Cluster by distance in X-Y //
    bool MinimumSeparation_xy(Float_t x_pos1, Float_t x_pos2, Float_t y_pos1, Float_t y_pos2);
    Double_t Separation_xy;
    Float_t x_step;
    Float_t xprev_step;
    Float_t y_step;
    Float_t yprev_step;


    
    // Temp arrays ( select, sort and cluster phases) //
    Int_t tempMax;
    TempArrays Temp;

};


#endif /* AnalysisSort_h */
