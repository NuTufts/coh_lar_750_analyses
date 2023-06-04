// c++ headers
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>
#include <stdlib.h>

// ROOT Headers
#include "TFile.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TChain.h"
#include "TApplication.h"
#include "TMath.h"

// CENNS Headers
#include "analysis/include/CENNSEvent.hh"
#include "analysis/include/CENNSChannel.hh"

int main( int nargs, char** argv ) {

  std::cout << "Plot Cherenkov Information From Hits" << std::endl;
  if ( nargs!=3 ) {
    std::cout << "usage: plot_cherenkov_hits [input file] [outputfile]" << std::endl;
    return 0;
  }
  
  std::string input_cenns_geant_file = argv[1];
  std::string output_root            = argv[2];

  std::ifstream f(output_root.c_str());

  if ( f.good() ) {
    std::cout << "Output file already exists. Remove first." << std::endl;
    return 0;
  }
  f.close();

  TApplication app( "app", &nargs, argv );
  
  // ------------------------------------------------------------
  // INPUT FILE
  
  TFile* finput = new TFile( input_cenns_geant_file.c_str(), "OPEN" );

  TTree* cenns = (TTree*)finput->Get("CENNS");
  CENNSEvent* eventdata = nullptr;
  TBranch* b_event = 0;
  cenns->SetBranchAddress("Event",&eventdata,&b_event);
  
  size_t nentries = cenns->GetEntries();

  std::cout << "Number of entries in file: " << nentries << std::endl;
  
  std::cout << "Input file: " << input_cenns_geant_file << " (length = " << input_cenns_geant_file.size() << " )" << std::endl;
  input_cenns_geant_file.erase(0, 56);
  std::cout << "Input file: " << input_cenns_geant_file << " (length = " << input_cenns_geant_file.size() << " )" << std::endl;

  if (input_cenns_geant_file.size() == 15) {
        input_cenns_geant_file.erase(input_cenns_geant_file.begin()+1, input_cenns_geant_file.end()-7);
        input_cenns_geant_file.erase(input_cenns_geant_file.begin()+3, input_cenns_geant_file.end());
  } else if (input_cenns_geant_file.size() == 16) {
        input_cenns_geant_file.erase(input_cenns_geant_file.begin()+2, input_cenns_geant_file.end()-7);
        input_cenns_geant_file.erase(input_cenns_geant_file.begin()+4, input_cenns_geant_file.end());
  }
  std::cout << "Energy: " << input_cenns_geant_file << std::endl;


  // ------------------------------------------------------------
  // OUTPUT FILE
  TFile* fout = new TFile( output_root.c_str(), "NEW" );

  // arrays: [detector][photon type]
  typedef enum { kScint=0, kCeren, kNumPhotonTypes } PhotonType_t;
  typedef enum { kPMT=0, kSiPM, kNumDetTypes } DetType_t;
  std::vector<std::string> photontype_v = {"scint", "ceren" };
  std::vector<std::string> dettype_v    = {"pmt"  , "sipm" };
  TRandom3* rand = new TRandom3();
  double randNum;
  // Primary colors are scintillation
  // Pastelle colors are for Cherenkov
  Color_t colors[2][2] = { {kBlue, kCyan},      // PMTs
                           {kRed,  kMagenta} }; // SiPMs

  // Detector Geometry
  float Radius     = 14.50*25.4; // mm
  float HalfHeight = 19.5*25.4;  // mm

  float histH = HalfHeight + 2*Radius+10;
  float histW = 2*Radius*3.14159+10;

  // Histogram for total number of hits
  TH1D* htotals[2][2];
  TH1D* hdet[2];
  for ( int idet_type=0; idet_type<kNumDetTypes; idet_type++ ) {
    float binmax = ( idet_type==kPMT ) ? 50.0*10e3 : 50.0*5;
    
    for ( int iphoton_type=0; iphoton_type<kNumPhotonTypes; iphoton_type++ ) {
      char histname[50];
      sprintf( histname, "htotals_%s_%s", dettype_v[idet_type].c_str(), photontype_v[iphoton_type].c_str() );
      htotals[idet_type][iphoton_type] = new TH1D( histname, ";num pe;counts",100, 0, binmax );
      htotals[idet_type][iphoton_type]->SetLineColor( colors[idet_type][iphoton_type] );
    }

    char histname[50];
    sprintf( histname, "htotals_%s", dettype_v[idet_type].c_str() );
    hdet[idet_type] = new TH1D( histname, ";num pe;counts",100, 0, binmax );
    hdet[idet_type]->SetLineColor( colors[idet_type][kScint] );
  }

  TTree* tottree = new TTree("hittotals","Hit Totals");
  float totals[2][2] = { {0,0}, {0,0} };

  tottree->Branch( "totals", totals, "totals[2][2]/F" );
  
  // make vectors to store hit totals for every event
  std::vector<int> pmt_scint_hit_totals;
  std::vector<int> sipm_ceren_hit_totals;
  
  std::string filename = "/cluster/tufts/wongjiradlabnu/bharri06/data/tot_data/TrueV_hittotals_data_10x/data_";
  filename += input_cenns_geant_file;
  filename += ".txt";
  std::cout << "Filename: " << filename << std::endl;
  
  // write yield files to store data in
  std::ofstream outputfile;
  outputfile.open(filename);

  if (input_cenns_geant_file.size() == 3) {
        input_cenns_geant_file.erase(input_cenns_geant_file.begin()+1, input_cenns_geant_file.end());
  } else if (input_cenns_geant_file.size() == 4) {
        input_cenns_geant_file.erase(input_cenns_geant_file.begin()+2, input_cenns_geant_file.end());
  } 

  std::cout << "Energy: " << input_cenns_geant_file << std::endl;
 
  /****** PLOT CHERENKOV POSITIONS *********/
  std::vector<double> energy;
  std::vector<double> true_z;
  std::vector<double> hit_totals[2];
  
  
  //Event Loop
  //
  //
  //
  for ( size_t ientry=0; ientry<nentries; ientry++ ) {

    std::cout << "[ENTRY " << ientry << "]" << std::endl;
    
    ULong_t bytes = cenns->GetEntry(ientry);
    if ( bytes==0 ) {
      // empty event
      continue;
    }

    std::vector<int> tb(2,0);

    for (int i=0; i<2; i++ )
      for (int j=0; j<2; j++ )
        totals[i][j] = 0.0;

    for ( auto const& chdata : eventdata->fvChannels ) {
      int channelid = chdata.GetChannelNum();

      int numhits = chdata.size();
      //std::cout << "  Channel[" << channelid << "] numhits=" << numhits << std::endl;
      for ( auto const& photonhit : chdata ) {

	// QE required variables
        // float eV = photonhit.energy_ev;
        // float nm = 1240/eV; // 1240 = nm * eV

        int is_ceren = photonhit.ischerenkov;
	
	//get z positon of hit
	double HitPosZ = photonhit.pos[2];
        
        if ( !is_ceren ) {
          if ( channelid<300 ) {
	    totals[kPMT][kScint]+=1;
	    if ( HitPosZ>0 ) { 
              tb[0]++;
	    } else { 
              tb[1]++;
            }
	  } else {
            totals[kSiPM][kScint]+=1;
	  }
	}
        else {
          if ( channelid<300 )
            totals[kPMT][kCeren]+=1;
          else
            totals[kSiPM][kCeren]+=1;
        }

      }// end of hit loop

    }//end of channel loop
    
    std::cout << "Out of channel loop." << std::endl;


    //keep track of total scintillation and cherenkov hits recieved for every event;
    pmt_scint_hit_totals.push_back( totals[kPMT][kScint]);
    sipm_ceren_hit_totals.push_back( totals[kSiPM][kCeren]);
	
    std::cout << "made it" << std::endl; 
    // fill histograms
    for (int i=0; i<2; i++ ) {
      for (int j=0; j<2; j++) {
        htotals[i][j]->Fill( totals[i][j] );
      }
    }
    hdet[kPMT]->Fill( totals[kPMT][kScint]+totals[kPMT][kCeren] );
    hdet[kSiPM]->Fill( totals[kSiPM][kScint]+totals[kSiPM][kCeren] );

    // fill tree
    tottree->Fill();

    energy.push_back( std::stod(input_cenns_geant_file) );
    true_z.push_back( eventdata->fGenZ );
    hit_totals[0].push_back( tb[0] );
    hit_totals[1].push_back( tb[1] );
     
  } // end of event loop
   
  std::cout << "Out of event loop." << std::endl;

  std::cout << "done with loop" << std::endl;
  
  for (int i=0; i<nentries; i++) {
  	outputfile << energy.at(i) << " " << true_z.at(i) << " " << hit_totals[0].at(i) << " " << hit_totals[1].at(i) << std::endl;
  }
  outputfile.close();

  fout->Write();
  fout->Close();
  finput->Close();

  
  return 0;
}
