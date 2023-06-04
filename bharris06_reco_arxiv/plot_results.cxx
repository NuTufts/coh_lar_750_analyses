// c++ headers
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
using namespace std;

// ROOT Headers
//#include "/usr/local/root/include/*.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TDialogCanvas.h"
#include "TInspectCanvas.h"
#include "TTree.h"
#include "TChain.h"
//#include "TApplication.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraph2DErrors.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TStyle.h"

int main( int nargs, char** argv ) {

  std::cout << "Correcting vertex using data written to text file by reco_vertex..." << std::endl;
  if ( nargs!=2 ) {
    std::cout << "usage: corr_vertex [outputfile]" << std::endl;
    return 0;
  }

  std::string output_root     = argv[1];

  std::ifstream f(output_root.c_str());

  if ( f.good() ) {
    std::cout << "Output file already exists. Remove first." << std::endl;
    return 0;
  }
  f.close();
  // ------------------------
  // CREATE OUTPUT ROOT FILES
  // ------------------------
  TFile* fout = new TFile( output_root.c_str(), "NEW" );
  // --------------------------------------------------
  // READ IN INITIAL VERTEX GUESS FROM INPUT TEXT FILES
  // --------------------------------------------------
  // create vectors to store data from input text files
  std::vector<double> tb[10][2];             //[top,bot]
  std::vector<double> true_vertex[10][3];   //[x,y,z]
  std::vector<double> reco_vertex[10][3];   //[x,y,z]
  std::vector<double> vertex_error[10];
  std::vector<double> true_dir[10][3];      //[x,y,z]
  std::vector<double> reco_dir[10][3];      //[x,y,z]
  std::vector<double> dir_error[10];  
  
  // loop through all input text files
  for ( int x=5; x<51; x+=5) {
    for ( int y=0; y<50; y++) {
      // open text files
      std::string line;
      std::string filename = "/cluster/tufts/wongjiradlabnu/bharri06/data/reco_data/final/data_";
      filename += std::to_string(x);
      filename += "_";
      filename += std::to_string(y);
      filename += ".txt";
      std::ifstream myfile( filename);
      if ( myfile.is_open()) {
        while ( getline( myfile, line)) {
	  // read in variables from file
	  // first grab the energy and make the energy index
	  int ei =  std::stoi( line.substr(0, line.find(' ')) )/5-1;
	  line.erase(0, line.find(' ')+1);
	  tb[ei][0].push_back( std::stod( line.substr(0, line.find(' '))));
	  line.erase(0, line.find(' ')+1);
	  tb[ei][1].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  true_vertex[ei][0].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  true_vertex[ei][1].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  true_vertex[ei][2].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  reco_vertex[ei][0].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
          reco_vertex[ei][1].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
          reco_vertex[ei][2].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  true_dir[ei][0].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  true_dir[ei][1].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  true_dir[ei][2].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  reco_dir[ei][0].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  reco_dir[ei][1].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  reco_dir[ei][2].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  dir_error[ei].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  // find the reconstructed vertex error
	  vertex_error[ei].push_back( sqrt( pow(reco_vertex[ei][0].back()-true_vertex[ei][0].back(),2) +
				 		      pow(reco_vertex[ei][1].back()-true_vertex[ei][1].back(),2) +
				 		      pow(reco_vertex[ei][2].back()-true_vertex[ei][2].back(),2) ) );
	  //remove the RARE bad event (1 in 50k)
	  //files seem to be getting more corrupted each time?? 
	  //need to make new events but done have time right now :(
	  if (isnan(dir_error[ei].back())) {
	    tb[ei][0].pop_back();
	    tb[ei][1].pop_back();
            true_vertex[ei][0].pop_back();
            true_vertex[ei][1].pop_back();
            true_vertex[ei][2].pop_back();
            reco_vertex[ei][0].pop_back();
            reco_vertex[ei][1].pop_back();
            reco_vertex[ei][2].pop_back();
	    true_dir[ei][0].pop_back();
	    true_dir[ei][1].pop_back();
	    true_dir[ei][2].pop_back();
	    reco_dir[ei][0].pop_back();
	    reco_dir[ei][1].pop_back();
	    reco_dir[ei][2].pop_back();
	    dir_error[ei].pop_back();
	    vertex_error[ei].pop_back();
	  }
	}  
        myfile.close();	
      } else { 
	std::cout << "Unable to open file: " << filename << std::endl;
      }
    }
  }
  

  // ---------------------------------------------
  // Print number of events per energy level
  // ---------------------------------------------
  for (int i=0; i<10; i++) {
    std::cout << i*5+5 << " MeV has " << reco_vertex[i][0].size() << " events." << std::endl;
  } 
 

 // --------------------------------------------------------
 // Find the mean and absolute vertex error at each energy 
 // --------------------------------------------------------
 std::vector<double> absolute_vertex_error(10,0);
 std::vector<double> mean_vertex_error(10,0);
 //make a list of values used in mean vertex error to make sd
 std::vector<double> list_vertex_error[10];
 // loop through every energy
 for (int e=0;e<10;e++) {
   // loop through every event at that energy
   int size=reco_vertex[e][2].size();
   for (int i=0;i<size;i++) {
     list_vertex_error[e].push_back(vertex_error[e][i]);
     mean_vertex_error[e]+=list_vertex_error[e].back();
   } //end of event loop 
   mean_vertex_error[e]=mean_vertex_error[e]/size;
   //find the max error of the best 68% of results
   std::sort(list_vertex_error[e].begin(), list_vertex_error[e].end());
   absolute_vertex_error[e]=list_vertex_error[e].at(0.68*size);
 } //end of energy loop

 // ----------------------------------------------------
 // CALCULATE THE STANDARD DEVIATION OF THE MEAN VERTEX ERROR 
 // ----------------------------------------------------
 std::vector<double> sd_mean_vertex_error(10,0);
 for (int e=0;e<10;e++) {
   //loop through every hit in that zone and sum the difference squared
   int size=reco_vertex[e][2].size();
   for (int i=0;i<size;i++) {
     sd_mean_vertex_error[e]+=pow(list_vertex_error[e].at(i)-mean_vertex_error[e],2);
   }
   // divide by total number of entries for that energy and take the sqrt
   sd_mean_vertex_error[e]=sd_mean_vertex_error[e]/size;
   sd_mean_vertex_error[e]=sqrt(sd_mean_vertex_error[e]);
 } // end of energy loop

 // --------------------------------------------------------
 // Find the mean and absolute direction error at each energy 
 // --------------------------------------------------------
 std::vector<double> absolute_dir_error(10,0);
 std::vector<double> mean_dir_error(10,0);
 std::vector<double> list_dir_error[10];
 //loop through every energy
 for (int e=0;e<10;e++) { 
   // loop through every event at that energy
   int size=reco_vertex[e][2].size();
   for (int i=0;i<size;i++) {
      list_dir_error[e].push_back(dir_error[e][i]);  
      mean_dir_error[e]+=list_dir_error[e].back();
   } //end of event loop
   mean_dir_error[e]=mean_dir_error[e]/size;
   //find the max error of the best 68% of results
   std::sort(list_dir_error[e].begin(), list_dir_error[e].end());
   absolute_dir_error[e]=list_dir_error[e].at(0.68*size);
 }
 
 // ----------------------------------------------------
 // CALCULATE THE STANDARD DEVIATION OF THE MEAN DIR ERROR 
 // ----------------------------------------------------
 std::vector<double> sd_mean_dir_error(10,0);
 for (int e=0;e<10;e++) {
   //loop through every hit in that zone and sum the difference squared
   int size=reco_vertex[e][2].size();
   for (int i=0;i<size;i++) {
     sd_mean_dir_error[e]+=pow(list_dir_error[e].at(i)-mean_dir_error[e],2);
   }
   // divide by total number of entries for that energy and take the sqrt
   sd_mean_dir_error[e]=sd_mean_dir_error[e]/size;
   sd_mean_dir_error[e]=sqrt(sd_mean_dir_error[e]);
 } // end of energy loop

  // --------------------------------------------------------------
  // Plot the mean and absolute vertex and direction error as a function o energy
  // -------------------------------------------------------------
  TCanvas *c2 = new TCanvas("c2","reco vertex error", 5000, 5000);
  c2->cd(); 
  c2->SetGrid();

  //make energies to plot
  std::vector<double> energies = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50};
  std::vector<double> blank_sd(10,0);
 
  std::cout << "Energy        Abs. spatial res.      Mean spatial res.         Abs. angular res.          Mean angular res." << std::endl; 
  for (int i=0;i<10;i++) {
    std::cout << energies[i] << " " << absolute_vertex_error[i] << " " << mean_vertex_error[i] << " " << absolute_dir_error[i] << " " << mean_dir_error[i] << std::endl;
  }

  TGraphErrors *mean_spatial_res = new TGraphErrors(10,energies.data(),mean_vertex_error.data(),blank_sd.data(),sd_mean_vertex_error.data());
  mean_spatial_res->SetTitle("Mean Spatial Resolution; Energy of e- Produced (MeV); Mean Error of Reconstructed Vertex (mm)");
  mean_spatial_res->SetMarkerColor(4);
  mean_spatial_res->SetMarkerStyle(21);
  mean_spatial_res->SetName("mean_spatial_res");
  mean_spatial_res->Write();

  TGraph *abs_spatial_res = new TGraph(10,energies.data(),absolute_vertex_error.data());
  abs_spatial_res->SetTitle("Absolute Spatial Resolution; Energy of e- Produced (MeV); Absolute Error of Reconstructed Vertex (mm)");
  abs_spatial_res->SetMarkerColor(4);
  abs_spatial_res->SetMarkerStyle(21);
  abs_spatial_res->SetName("abs_spatial_res");
  abs_spatial_res->Write();

  TGraphErrors *mean_angular_res = new TGraphErrors(10,energies.data(),mean_dir_error.data(),blank_sd.data(),sd_mean_dir_error.data());
  mean_angular_res->SetTitle("Mean Angular Resolution; Energy of e- Produced (MeV); Mean Error of Reconstructed Direction (deg)");
  mean_angular_res->SetMarkerColor(2);
  mean_angular_res->SetMarkerStyle(21);
  mean_angular_res->SetName("mean_anular_res");
  mean_angular_res->Write();  
  
  TGraph *abs_angular_res = new TGraph(10,energies.data(),absolute_dir_error.data());
  abs_angular_res->SetTitle("Absolute Angular Resolution; Energy of e- Produced (MeV); Absolute Error of Reconstructed Direction (deg)");
  abs_angular_res->SetMarkerColor(2);
  abs_angular_res->SetMarkerStyle(21);
  abs_angular_res->SetName("abs_angular_res");
  abs_angular_res->Write();

  //------------------------------------------------------------------------------------------------------------
  //Find correlation betwen total number of scintillation hits and the true energy of the e- created 
  //------------------------------------------------------------------------------------------------------------								
  
  //plot all results to see
  std::stringstream energyreco_name;
  energyreco_name << "energyreco";
  TH2D energyreco( energyreco_name.str().c_str(), "for correlation purposes; Total number of scintillation hits observed; True energy of e-", 50000, 0, 50000, 100, 0, 55);

  for (int ei=0;ei<10;ei++) {
    int size = tb[ei][0].size();
    for (int i=0;i<size;i++) {
      energyreco.Fill(tb[ei][0][i]+tb[ei][1][i],5*ei+5);
    }
  }
  energyreco.Write();
  
  //----------------------------------------------------------------------------
  //plot the mean total scintillation hits for each energy and fit a line to it
  //--------------------------------------------------------------------------
  std::vector<double> mean_tsh(10,0);
  for (int ei=0;ei<10;ei++) {
    int size=tb[ei][0].size();
    for (int i=0;i<size;i++) {
      mean_tsh[ei]+=tb[ei][0][i]+tb[ei][1][i];
    }
    mean_tsh[ei]=mean_tsh[ei]/size;
    //std::cout << mean_tsh[ei] << " " << std::endl;
  }

  //----------------------------
  //try median instead of mean
  //----------------------------
  std::vector<double> list_tsh[10];
  std::vector<double> median_tsh(10,0);
  for (int ei;ei<10;ei++) {
    int size=tb[ei][0].size();
    for (int i=0;i<size;i++) {
      list_tsh[ei].push_back(tb[ei][0][i]+tb[ei][1][i]);
    }
    std::sort(list_tsh[ei].begin(), list_tsh[ei].end());
    median_tsh[ei]=list_tsh[ei][0.5*(size-1)];
  }


  TGraph *mtsh = new TGraph(10,median_tsh.data(),energies.data());
  mtsh->SetTitle("for correlation purposes; Median number of scintillation hits observered; True energy of e-"); 
  mtsh->SetName("med_tsh");
  mtsh->SetMarkerColor(3);
  mtsh->SetMarkerStyle(21);

 // Fit Code
    double xmin = median_tsh[0];
    double xmax = median_tsh[median_tsh.size()-1];
    TF1 myxFunc("Polynomial3rdDegree","pol3",xmin,xmax);
    mtsh->Fit(&myxFunc);
    double x0 = myxFunc.GetParameters()[0];
    double x1 = myxFunc.GetParameters()[1];
    double x2 = myxFunc.GetParameters()[2];
    double x3 = myxFunc.GetParameters()[3];
    //double x4 = myxFunc.GetParameters()[4];
    //double x5  = myxFunc_fifth.GetParameters()[5];
    TGraph xline(0);
    int xi =0;
    for (int x_coord = xmin; x_coord <= xmax; x_coord++){
      double y= x0*pow(x_coord,3) + x1*pow(x_coord,2) + x2*x_coord + x3;
      xline.SetPoint(xi,x_coord,y);
      xi++;
    }
    mtsh->Write();

    std::cout << x3 << " " << x2 << " " << x1 << " " << x0 << std::endl;
    std::cout << 4.84643*(1/pow(10,13)) << " " << -4.8728*(1/pow(10,10)) << " " << 0.00202442 << " " << -0.0569654  << std::endl;
  

  //----------------------------------------------------------------------------------
  //See how well we can reconstruct the energy using the mean total scintillation hits
  //---------------------------------------------------------------------------------

  std::vector<double> abs_eerror(10,0);
  std::vector<double> mean_eerror(10,0);
  std::vector<double> list_eerror[10];
  //loop through every event and reco its energy, then keep track of the error
  for (int ei=0;ei<10;ei++) {
    int size = tb[ei][0].size();
    for (int i=0;i<size;i++) {
      double tot = tb[ei][0][i]+tb[ei][1][i];
      double reco_e = x3*pow(tot,3)+x2*pow(tot,2)+x1*tot+x0;
      ////beam we will use only generates up to 52.8 MeV
      //if (reco_e > 52.8) {
      //  reco_e=52.8;
      //}
      double true_e = 5*ei+5;
      double error = abs(reco_e - true_e);
      mean_eerror[ei]+=error/true_e;
      list_eerror[ei].push_back(error/true_e);
    }
    mean_eerror[ei]=mean_eerror[ei]/size;
    std::sort(list_eerror[ei].begin(), list_eerror[ei].end());
    abs_eerror[ei]=list_eerror[ei].at(0.68*size);
    std::cout << mean_eerror[ei] << " " <<  std::endl;
  }
  

  TGraph *abs_energy_res = new TGraph(10,energies.data(),abs_eerror.data());
  abs_energy_res->SetTitle("Absolute Energy Resolution; Energy of e- Produced (MeV); Absolute Percent Error of Reconstructed Direction (MeV)");
  abs_energy_res->SetMarkerColor(3);
  abs_energy_res->SetMarkerStyle(21);
  abs_energy_res->SetName("abs_energy_res");
  abs_energy_res->Write(); 
  
   // -----------------------------------------------------
 // Make 1D histo of the energy,  vertex and direction error at 5, 20 and 50 MeV
 // -----------------------------------------------------

 int ee=0;
 const int n =reco_vertex[ee][2].size();
 TH1D *verror5mev = new TH1D("verror5mev","blank8", n, 0, 300);
 TH1D *derror5mev = new TH1D("derror5mev","blank7", n, 0, 200);
 TH1D *eerror5mev = new TH1D("eerror5mev","blank9", n, 0, 1);
 TH1D *scint5mev = new TH1D("scint5mev","blank10", n, 0, 5000);
 for (int i=0; i<n; i++) {
  verror5mev->Fill(vertex_error[ee][i]);
  derror5mev->Fill(dir_error[ee][i]);
  eerror5mev->Fill(list_eerror[ee][i]);
  scint5mev->Fill(tb[ee][0][i]+tb[ee][1][i]);
 }

 ee=3;
 const int m = reco_vertex[ee][2].size();
 TH1D *verror20mev = new TH1D("verror20mev","blank88", m, 0, 300);
 TH1D *derror20mev = new TH1D("derror20mev","blank77", m, 0, 200);
 for (int i=0; i<m; i++) {
   verror20mev->Fill(vertex_error[ee][i]);
   derror20mev->Fill(dir_error[ee][i]);
 }

 ee=9;
 const int o = reco_vertex[ee][2].size();
 TH1D *verror50mev = new TH1D("verror50mev","blank888", o, 0, 300);
 TH1D *derror50mev = new TH1D("derror50mev","blank777", o, 0, 200);
 TH1D *eerror50mev = new TH1D("eerror50mev","blank999", o, 0, 1);
 TH1D *scint50mev = new TH1D("scint50mev","blank10", n, 0, 35000);
 for (int i=0; i<o; i++) {
   verror50mev->Fill(vertex_error[ee][i]);
   derror50mev->Fill(dir_error[ee][i]);
   eerror50mev->Fill(list_eerror[ee][i]);
   scint50mev->Fill(tb[ee][0][i]+tb[ee][1][i]);
 }

 verror5mev->SetTitle("Error of Reconstructed Vertex at 5 MeV; Reconstructed Vertex Error (mm); Number of Events");
 verror5mev->Write();
 verror20mev->SetTitle("Error of Reconstructed Vertex at 20 MeV; Reconstructed Vertex Error (mm); Number of Events");
 verror20mev->Write();
 verror50mev->SetTitle("Error of Reconstructed Vertex at 50 MeV; Reconstructed Vertex Error (mm); Number of Events");
 verror50mev->Write();

 derror5mev->SetTitle("Error of Reconstructed Direction at 5 MeV; Reconstructed Direction Error (deg); Number of Events");
 derror5mev->Write();
 derror20mev->SetTitle("Error of Reconstructed Direction at 20 MeV; Reconstructed Direction Error (deg); Number of Events");
 derror20mev->Write();
 derror50mev->SetTitle("Error of Reconstructed Direction at 50 MeV; Reconstructed Direction Error (deg); Number of Events");
 derror50mev->Write();
 
 eerror5mev->SetTitle("Percent Error of Reconstructed Energy at 5 MeV; Reconstructed Energy Percent Error (MeV); Number of Events");
 eerror5mev->Write();
 eerror50mev->SetTitle("Percent Error of Reconstructed Energy at 50 MeV; Reconstructed Energy Percent Error (MeV); Number of Events");
 eerror50mev->Write();

 scint5mev->SetTitle("Number of Scintillation Hits Observed at 5 Mev; Number of Scintillation Hits Observed; Number of Events");
 scint5mev->Write();
 scint50mev->SetTitle("Number of Scintillation Hits Observed at 50 Mev; Number of Scintillation Hits Observed; Number of Events");
 scint50mev->Write();
 
  






  //-------------close output root file--------------------

  fout->Close();

  return 0;
}
