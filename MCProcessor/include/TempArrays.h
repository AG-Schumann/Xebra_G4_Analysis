//
//  TempArrays.h
//  Clustering
//
//  Created by Yanina on 07.09.18.
//

#ifndef TempArrays_h
#define TempArrays_h

#include <stdio.h>
#include <vector>
#include <RtypesCore.h>
#include <stddef.h>

using namespace std;

struct TempArrays{
    Float_t Temp_x, Temp_y, Temp_z, Temp_etot;
    
    TempArrays() {}
    
    bool operator < ( const TempArrays& struct1) const {
      
       return  Temp_z < struct1.Temp_z;
       
    }

};



#endif /* TempArrays_h */
