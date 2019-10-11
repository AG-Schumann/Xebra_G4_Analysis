//
// MainProcessor.cc
// XeBRA MC Processor
//
// Created by Patricia on 27.11.18. 
// Modified for DBSCAN PreClustering by F. Kuger (Mai 2019)
//

// Processor Libraries
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include "TBranch.h"
#include <algorithm>
#include <vector>
#include <ctime>
#include "DBSCAN.hpp"
#include <Eigen/Eigen/Core>

//Processor Headers
#include "MainProcessor.h"
#include "TempArrays.h"


using namespace std;
using namespace clustering; 
using namespace Eigen; 

MainProcessor::MainProcessor(TString  InputFilename,TString  OutputFilename) : inTree(0){;}

MainProcessor::~MainProcessor(){;}

void MainProcessor::Loop(TString  InputFilename, TString  OutputFilename) {

/////////////////////////////////////////////////////////////////
///// Inititalize ROOT Input / Output 
/////////////////////////////////////////////////////////////////
	
	// Load Input file and Input Tree //
	TFile *inFile = new TFile(InputFilename, "READ");
	TTree *inTree = (TTree*)inFile->Get("events/events");
	
	// Load Output file and Output Tree //
	TFile *fOutTree = new TFile(OutputFilename,"recreate");
	TTree *outTree = new TTree("events","ProcessingFile");

	Long64_t nentries = inTree -> GetEntriesFast();
	
	inTree->SetBranchAddress("eventid", &eventid);
	inTree->SetBranchAddress("etot", &etot);
	inTree->SetBranchAddress("nsteps", &nsteps);
	inTree->SetBranchAddress("e_pri", &e_pri);
	inTree->SetBranchAddress("xp_pri", &xp_pri);
	inTree->SetBranchAddress("yp_pri", &yp_pri);
	inTree->SetBranchAddress("zp_pri", &zp_pri);
	inTree->SetBranchAddress("xp", &xp);
	inTree->SetBranchAddress("yp", &yp);
	inTree->SetBranchAddress("zp", &zp);
	inTree->SetBranchAddress("ed", &ed);
	inTree->SetBranchAddress("type_pri", &type_pri);
 	inTree->SetBranchAddress("type", &type);
   	inTree->SetBranchAddress("creproc", &creaproc);
	inTree->SetBranchAddress("edproc", &edproc);
	inTree->SetBranchAddress("LScint_etot", &LScint_etot);


	//output tree variables : set branches
	outTree->Branch("Xp", &Xp);
	outTree->Branch("Yp", &Yp);
	outTree->Branch("Zp", &Zp);
	outTree->Branch("Etot", &Etot);
	outTree->Branch("Xp_RMS", &Xp_RMS);
	outTree->Branch("Yp_RMS", &Yp_RMS);
	outTree->Branch("Zp_RMS", &Zp_RMS);
	outTree->Branch("type_primary", "vector<string>", &type_primary);
	outTree->Branch("type_particle", "vector<string>", &type_particle);
	outTree->Branch("process", "vector<string>", &process);
	outTree->Branch("nScat", &nScat, "nScat/I");
	outTree->Branch("xpri", &xpri, "xpri/F");
	outTree->Branch("ypri", &ypri, "ypri/F");
	outTree->Branch("zpri", &zpri, "zpri/F");
	outTree->Branch("epri", &epri, "epri/F");
	outTree->Branch("Edtot", &Edtot);
	outTree->Branch("LScintEdtot", &LScintEdtot);
   
/////////////////////////////////////////////////////////////////
///// Inititalize the DBSCAN PreClustering Algorithm   
/////////////////////////////////////////////////////////////////

   	// Parameters for DBSCAN Preclustering 
	double eps = 5.0;	//epsilon distance in mm
	int min_elems = 1;

   	// Cast an instance of the DBSCAN Template Class using the predetermined parameters for epsilon and Min_elems
	DBSCAN<Eigen::VectorXd, Eigen::MatrixXd> *myDBSCAN = new DBSCAN <Eigen::VectorXd, Eigen::MatrixXd>(eps, min_elems);	
	myDBSCAN->init(eps, min_elems);	



///////////////////////////////////////////////////////////////////
///// Starting Main Loop over tree entries 
/////////////////////////////////////////////////////////////////

	
	for(Long64_t jentry = 0; jentry<nentries; jentry++){

		inTree->GetEntry(jentry);

		if ( jentry % 1000 ==0 )  cout << jentry << " Entries  " << endl;

		//We apply the cuts before calling the AnalysisSort function
		if( nsteps == 0 ) continue;	  //remove events with no steps 
		if( etot == 0 ) continue;	  //remove events with no energy depositions
		

   		/////////////////////////////////////////////////////////////////
   		///// PreClustering with DBScan
   		/////////////////////////////////////////////////////////////////
		
		//Convert E_D into Input Matrix (+ Boundary points)
		int rows= xp->size() +2 ;
		int columns = 3;
		Eigen::MatrixXd myEDepositions(rows, columns); 
		for (int i= 0; i< (rows-2); i++){ myEDepositions(i,0) = (*xp)[i]; myEDepositions(i,1) = (*yp)[i]; myEDepositions(i,2) = (*zp)[i];}
		for (int k =0; k< 3; k++) {myEDepositions(rows-2,k) = 5000.; myEDepositions(rows-1,k) = -5000.;}
		
		//Run DBSCAN and retrieve PreCluster_ID - Label
		myDBSCAN->fit(myEDepositions);
		vector<int> myLabels = myDBSCAN->get_labels();
		myDBSCAN->reset(); 

		//Calculate PreCluster Information & Format Output Vectors

		int nPreCluster_temp = *(max_element(myLabels.data(), myLabels.data() + myLabels.size())) -1;
	
		vector<float> x_PreClusters_temp (nPreCluster_temp, 0.0);
		vector<float> y_PreClusters_temp (nPreCluster_temp, 0.0);
		vector<float> z_PreClusters_temp (nPreCluster_temp, 0.0);
		vector<float> e_PreClusters_temp (nPreCluster_temp, 0.0);

		vector<float> x_RMS_PreClusters_temp (nPreCluster_temp, 0.0);
		vector<float> y_RMS_PreClusters_temp (nPreCluster_temp, 0.0);
		vector<float> z_RMS_PreClusters_temp (nPreCluster_temp, 0.0);

		// Calculate weighted positions and total energy
		for (int i= 0; i< (rows-2); i++){
		int clusterID = myLabels[i];
			if((*ed)[i] == 0 ){(*ed)[i] = 0.001;}
		x_PreClusters_temp[clusterID] = x_PreClusters_temp[clusterID] + (*xp)[i] * (*ed)[i];
		y_PreClusters_temp[clusterID] = y_PreClusters_temp[clusterID] + (*yp)[i] * (*ed)[i];
		z_PreClusters_temp[clusterID] = z_PreClusters_temp[clusterID] + (*zp)[i] * (*ed)[i];
		e_PreClusters_temp[clusterID] = e_PreClusters_temp[clusterID] + (*ed)[i];
		}
		for (int j= 0; j< nPreCluster_temp; j++){
		x_PreClusters_temp[j] = x_PreClusters_temp[j] / e_PreClusters_temp[j];
		y_PreClusters_temp[j] = y_PreClusters_temp[j] / e_PreClusters_temp[j];		
		z_PreClusters_temp[j] = z_PreClusters_temp[j] / e_PreClusters_temp[j];
		}

		// Calculate the RMS
		for (int i= 0; i< (rows-2); i++){
		int clusterID = myLabels[i];
		x_RMS_PreClusters_temp[clusterID] = x_RMS_PreClusters_temp[clusterID] + (pow((*xp)[i] - x_PreClusters_temp[clusterID],2) * (*ed)[i]);
		y_RMS_PreClusters_temp[clusterID] = y_RMS_PreClusters_temp[clusterID] + (pow((*yp)[i] - y_PreClusters_temp[clusterID],2) * (*ed)[i]);
		z_RMS_PreClusters_temp[clusterID] = z_RMS_PreClusters_temp[clusterID] + (pow((*zp)[i] - z_PreClusters_temp[clusterID],2) * (*ed)[i]);
		}
		for (int j= 0; j< nPreCluster_temp; j++){
		x_RMS_PreClusters_temp[j] = sqrt(x_RMS_PreClusters_temp[j] / e_PreClusters_temp[j]);
		y_RMS_PreClusters_temp[j] = sqrt(y_RMS_PreClusters_temp[j] / e_PreClusters_temp[j]) ;		
		z_RMS_PreClusters_temp[j] = sqrt(z_RMS_PreClusters_temp[j] / e_PreClusters_temp[j]) ;
		}

		

   		/////////////////////////////////////////////////////////////////
   		///// Removal of PreClusters with E < 0.1 keV
   		/////////////////////////////////////////////////////////////////
		double E_threshold = 0.1;
		int n_remove = 0;
		for (int i = 0; i<nPreCluster_temp; i++){
			if(e_PreClusters_temp[i] <= E_threshold) {n_remove = n_remove + 1;}
		}
				
		int nPreCluster = nPreCluster_temp - n_remove;
			
		vector<float> x_PreClusters (nPreCluster, 0.0);
		vector<float> y_PreClusters (nPreCluster, 0.0);
		vector<float> z_PreClusters (nPreCluster, 0.0);
		vector<float> e_PreClusters (nPreCluster, 0.0); 

		vector<float> x_RMS_PreClusters (nPreCluster, 0.0);
		vector<float> y_RMS_PreClusters (nPreCluster, 0.0);
		vector<float> z_RMS_PreClusters (nPreCluster, 0.0);
			
		n_remove = 0;
		for (int j= 0; j< nPreCluster_temp; j++){
			if(e_PreClusters_temp[j] <= E_threshold){n_remove = n_remove +1; }
			if(e_PreClusters_temp[j] > E_threshold){
			x_PreClusters[j - n_remove] = x_PreClusters_temp[j];
			y_PreClusters[j - n_remove] = y_PreClusters_temp[j];
			z_PreClusters[j - n_remove] = z_PreClusters_temp[j];
			e_PreClusters[j - n_remove] = e_PreClusters_temp[j];

			x_RMS_PreClusters[j - n_remove] = x_RMS_PreClusters_temp[j];
			y_RMS_PreClusters[j - n_remove] = y_RMS_PreClusters_temp[j];
			z_RMS_PreClusters[j - n_remove] = z_RMS_PreClusters_temp[j];
			}
		}

/*
		//Just some Terminal PrintOuts for debugging
		if ( jentry % 1000 ==235 )  {
		cout << jentry << " Entries  " << endl;

		std::cout << "Here is the matrix:\n" << myEDepositions << std::endl;
	
		std::cout << "Here is the Label Vector:\n" ;
		for (auto i = myLabels.begin(); i != myLabels.end(); ++i) {std::cout << *i << ' ';}
		std::cout << endl;

		std::cout << "number of preclusters:  " << nPreCluster << std::endl;
		for (int i = 0; i < nPreCluster; i++){ std::cout << x_PreClusters[i] << ' ' << y_PreClusters[i] << ' '<< z_PreClusters[i] << ' '<< e_PreClusters[i] << ' '<< endl;}
		}
*/	

   		/////////////////////////////////////////////////////////////////
   		///// Quanta generation & Propagation (NEST) --> To be done
   		/////////////////////////////////////////////////////////////////


   		/////////////////////////////////////////////////////////////////
   		///// Sensor Response --> To be done
   		/////////////////////////////////////////////////////////////////



   		/////////////////////////////////////////////////////////////////
   		///// Define Root return
   		/////////////////////////////////////////////////////////////////

		Xp = x_PreClusters;
		Yp = y_PreClusters;
		Zp = z_PreClusters;
		Etot = e_PreClusters;

		Xp_RMS = x_RMS_PreClusters;
		Yp_RMS = y_RMS_PreClusters;
		Zp_RMS = z_RMS_PreClusters;

		nScat = nPreCluster;
		type_primary = type_pri;
		type_particle = type;
		xpri = xp_pri;
		ypri = yp_pri;
		zpri = zp_pri;
		epri = e_pri;
		process = creaproc;
		Edtot = etot;
		LScintEdtot = LScint_etot;
		//// use if condition for single scatter cut:
		//if(nScat == 1){
		outTree -> Fill();
		//}

	
	} //end of the loop over the number of entries in the tree 



	fDebug = 0;
	if (fDebug == 1 )	outTree -> Print();
	
	// Version and Timestamp (WARNING, LOCAL TIME, NO UTC) //
	time_t now = time(0);char* dt = ctime(&now);TNamed time("time",dt);TNamed version("version","MC processor v1.0");
	time.Write();version.Write();
	
	outTree -> Write();
	fOutTree -> Close();

}



