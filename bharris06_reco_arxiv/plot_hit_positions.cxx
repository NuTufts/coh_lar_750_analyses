// c++ headers
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <sstream>

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
#include "/cluster/tufts/wongjiradlab/bharri06/cenns10geant4/analysis/include/CENNSEvent.hh"
#include "/cluster/tufts/wongjiradlab/bharri06/cenns10geant4/analysis/include/CENNSChannel.hh"

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
  
  // ------------------------------------------------------------
  // OUTPUT FILE
  TFile* fout = new TFile( output_root.c_str(), "NEW" );

  std::cout << "output file made" << std::endl;

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

  /****** PLOT CHERENKOV POSITIONS *********/
  //Event Loop
  for ( size_t ientry=0; ientry<20; ientry++ ) {

    std::cout << "[ENTRY " << ientry << "]" << std::endl;
    
    ULong_t bytes = cenns->GetEntry(ientry);
    if ( bytes==0 ) {
      // empty event
      continue;
    }

    std::vector<double> hit_x[2];
    std::vector<double> hit_y[2];
    std::vector<double> hit_z[2];

    std::stringstream hposz_name;
    hposz_name << "hposz_entry" << ientry;
    TH1D hposz( hposz_name.str().c_str(), ";z", 2000, -1000, 1000 );

    std::stringstream hposxy_top_name;
    std::stringstream hposxy_bot_name;    
    hposxy_top_name << "hposxy_top_entry" << ientry;
    hposxy_bot_name << "hposxy_bot_entry" << ientry;    

    // Scintillation Hits
    TH2D hposxy_top( hposxy_top_name.str().c_str(), ";x;y", 100, -400, 400, 100, -400, 400 );
    TH2D hposxy_bot( hposxy_bot_name.str().c_str(), ";x;y", 100, -400, 400, 100, -400, 400 );

    std::stringstream detview_scint_name;
    std::stringstream detview_ceren_name;
    detview_scint_name << "detview_scint_entry" << ientry;
    detview_ceren_name << "detview_ceren_entry" << ientry;
    
    TH2D detview_scint( detview_scint_name.str().c_str(), "", 500, -histW/2, histW/2, 500, -histH, histH);
    TH2D detview_ceren( detview_ceren_name.str().c_str(), "", 500, -histW/2, histW/2, 500, -histH, histH);
   


    std::vector<double> sorted_hit_time;

    // channel loop
    for ( auto const& chdata : eventdata->fvChannels ) {
      int channelid = chdata.GetChannelNum();
      int numhits = chdata.size();
      // hit loop
      for ( auto const& photonhit : chdata ) {
            // Get time of hit
            double HitTime = photonhit.time;
            int is_ceren = photonhit.ischerenkov;
            if (!is_ceren) {
                sorted_hit_time.push_back( HitTime );
            }
      }// end of hit loop
    }//end of channel loop

    std::sort( sorted_hit_time.begin(), sorted_hit_time.end());
    double first_hit_time = sorted_hit_time[0];


    for ( auto const& chdata : eventdata->fvChannels ) {
      int channelid = chdata.GetChannelNum();

      int numhits = chdata.size();
      std::cout << "  Channel[" << channelid << "] numhits=" << numhits << std::endl;
      for ( auto const& photonhit : chdata ) {

	// QE required variables
        float eV = photonhit.energy_ev;
        float nm = 1240/eV; // 1240 = nm * eV

        int is_ceren = photonhit.ischerenkov;
        
	double HitTime = photonhit.time;

	// Get hit positions
	double HitPosX = photonhit.pos[0];
	double HitPosY = photonhit.pos[1];
	double HitPosZ = photonhit.pos[2];
        double HitPosR = sqrt(HitPosX*HitPosX+HitPosY*HitPosY);

        if ( !is_ceren ) {
	  if (HitTime < first_hit_time+25) {
          	hit_x[0].push_back( HitPosX );
          	hit_y[0].push_back( HitPosY );
          	hit_z[0].push_back( HitPosZ );
          }
	}
        else {
          hit_x[1].push_back( HitPosX );
          hit_y[1].push_back( HitPosY );
          hit_z[1].push_back( HitPosZ );
        }

        hposz.Fill( HitPosZ );
      
	//if (HitTime < first_hit_time+25)
        if (channelid<114)
          if ( HitPosZ>0)
            hposxy_top.Fill( HitPosX, HitPosY );
          else
            hposxy_bot.Fill( HitPosX, HitPosY );
        // detector view
        float x,y;

        if (HitPosR < Radius-1.0 ) {
          // end caps
          float rot_theta = ( HitPosZ>0 ) ? -TMath::Pi()/2 : TMath::Pi()/2;
          x = HitPosX*cos(-rot_theta)  + HitPosY*sin(-rot_theta);
          y = -HitPosX*sin(-rot_theta) + HitPosY*cos(-rot_theta);
          if ( HitPosZ>0 )
            y += HalfHeight + Radius;
          else
            y -= (HalfHeight+Radius);
          
        }
        else {
          // wall
          x = 2*Radius*3.14159*(atan2( HitPosY, HitPosX )/(2*3.14159));
          y = HitPosZ;
        }
        if ( !is_ceren )
          detview_scint.Fill( x, y );
        else
	  // only take sipm ceren hits
	  if (channelid>114) {
            detview_ceren.Fill( x, y );
 	  }
      }// end of hit loop

    }//end of channel loop

    std::cout << "make 3d plot" << std::endl;
    if ( hit_x[0].size()>0 ) {
      TGraph2D* hit3d_scint = new TGraph2D( hit_x[0].size(), hit_x[0].data(), hit_y[0].data(), hit_z[0].data() );
      std::stringstream graphname_scint;
      graphname_scint << "hit3d_scint_entry" << ientry;
      hit3d_scint->Write( graphname_scint.str().c_str() );
      delete hit3d_scint;
    }

    if ( hit_x[1].size()>0 ) {
      TGraph2D hit3d_ceren( hit_x[1].size(), hit_x[1].data(), hit_y[1].data(), hit_z[1].data() );            
      std::stringstream graphname_ceren;
      graphname_ceren << "hit3d_ceren_entry" << ientry;
      hit3d_ceren.Write( graphname_ceren.str().c_str() );
    }

    hposz.Write();

    // make axis before writing plot to output root file
    hposxy_top.SetTitle("Scintillation Hits on Top End Cap; X (mm); Y (mm); Number of Hits");
    hposxy_top.Write();
    hposxy_bot.SetTitle("Scintillation Hits on Bottom End Cap; X (mm); Y (mm); Number of Hits");
    hposxy_bot.Write();

    detview_scint.Write();
    detview_ceren.SetTitle("Cherenkov Hits; X (mm); Y (mm); Number of Hits");
    detview_ceren.Write();

  } // end of event loop

  std::cout << "done with loop" << std::endl;
  
  //app.Run();

  fout->Close();
  finput->Close();


 

  
  return 0;
}
