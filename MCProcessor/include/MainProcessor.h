//
// MainProcessor.h
// XeBRA MC Processor
//
// Created by Patricia on 27.11.18.
//

#ifndef MainProcessor_h
#define MainProcessor_h

#include <stdio.h>
#include <vector>
#include <TTree.h>
#include <stddef.h>
#include <string>
#include "TempArrays.h"
#include "AnalysisSort.h"

using namespace std;

class MainProcessor{
public:
	// Constructor & Destructor //
	MainProcessor(TString  InputFilename,TString  OutputFilename);
	~MainProcessor();
	
	// Debug //
	Int_t fDebug;

	Long64_t nbytes = 0, nb = 0;

	// Input Tree //
	static TString InputFilename;
	TTree *inTree;
	Int_t		 eventid;
	Float_t		 etot;
	Int_t		 nsteps;
	Float_t		 e_pri;
	Float_t		 xp_pri;
	Float_t		 yp_pri;
	Float_t		 zp_pri;
	vector<string>	*type_pri=0;
	vector<string>	*type=0;
	vector<string>	*creaproc=0;
	vector<string>  *edproc=0;
	vector<float>   *xp = 0;
	vector<float>   *yp = 0;
	vector<float>   *zp = 0;
	vector<float>   *ed = 0;
	Float_t		LScint_etot;		

	// Output Tree //
	TFile *fOutTree;
	static TString OutputFilename;
	TTree *outTree;
	vector<string>	*type_primary;
	vector<string>	*type_particle;
	vector<string>	*process;
	vector<float>   Xp;
	vector<float>   Yp;
	vector<float>   Zp;
	vector<float>   Etot;
	vector<float>   Xp_RMS;
	vector<float>   Yp_RMS;
	vector<float>   Zp_RMS;
	Int_t		nScat;
	Float_t		xpri;
	Float_t		ypri;
	Float_t		zpri;
	Float_t		epri;
	Float_t		Edtot;
	Float_t		LScintEdtot;

	// Loop over events //
	void Loop(TString InputFilename, TString  OutputFilename);
	
};



#endif /* MainProcessor_h */
