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
  std::vector<double> true_vertex[10][4];   //[x,y,z,r]
  std::vector<double> reco_vertex[10][4];   //[x,y,z,r]
  std::vector<double> error[10][2];          //[r,z]
  
  //variables that keep track of number of events
  std::vector<std::vector<int>> events_per_file(10, std::vector<int> (50,0));
  std::vector<int> events_excluded(10,0);
  // loop through all input text files
  for ( int x=5; x<51; x+=5) {
    for ( int y=0; y<50; y++) {
      // open text files
      std::string line;
      std::string filename = "/cluster/tufts/wongjiradlabnu/bharri06/data/combined_reco/reco103/data_";
      filename += std::to_string(x);
      filename += "_";
      filename += std::to_string(y);
      filename += ".txt";
      std::ifstream myfile( filename);
      if ( myfile.is_open()) {
        while ( getline( myfile, line)) {
	  //std::cout << "Line: " << line << std::endl;
	  // read in variables from file
	  int energy_index =  std::stoi( line.substr(0, line.find(' ')) )/5-1;
	  line.erase(0, line.find(' ')+1);
	  tb[energy_index][0].push_back( std::stod( line.substr(0, line.find(' '))));
	  line.erase(0, line.find(' ')+1);
	  tb[energy_index][1].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  true_vertex[energy_index][0].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  true_vertex[energy_index][1].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  true_vertex[energy_index][2].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
	  reco_vertex[energy_index][0].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
          reco_vertex[energy_index][1].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);
          reco_vertex[energy_index][2].push_back( std::stod( line.substr(0, line.find(' '))));
          line.erase(0, line.find(' ')+1);	
	  //cout data grabbed from line
          //std::cout << "Cout: " << energy_index << " " << tb[energy_index][0].back() << " " << tb[energy_index][1].back() << " " << true_vertex[energy_index][0].back() << " " << true_vertex[energy_index][1].back() << " " << true_vertex[energy_index][2].back() << " " << reco_vertex[energy_index][0].back() << " " << reco_vertex[energy_index][1].back() << " " << reco_vertex[energy_index][2].back() << std::endl;
	  // find the z and r error of the reconstructed vertex
	  true_vertex[energy_index][3].push_back( sqrt( pow(true_vertex[energy_index][0].back(),2) + pow(true_vertex[energy_index][1].back(),2) ) );
	  reco_vertex[energy_index][3].push_back( sqrt( pow(reco_vertex[energy_index][0].back(),2) + pow(reco_vertex[energy_index][1].back(),2) ) );
	  error[energy_index][0].push_back( reco_vertex[energy_index][3].back() - true_vertex[energy_index][3].back() );
	  error[energy_index][1].push_back( reco_vertex[energy_index][2].back() - true_vertex[energy_index][2].back() );
	  //std::cout << "reco z: " << reco_vertex[energy_index][2].back() << " true z: " << true_vertex[energy_index][2].back() << std::endl;
	  //std::cout << "z error: " << error[energy_index][1].back() << std::endl;
	  
	  //if statement below is not active right now
	  // if the reconstructed radius is not in the middle of the detector, dont consider it
	  if (reco_vertex[energy_index][2].back() > 200 && reco_vertex[energy_index][2].back() < -200) {
            //get rid of the event you just added
	    tb[energy_index][0].pop_back();
	    tb[energy_index][1].pop_back();
	    true_vertex[energy_index][0].pop_back();
	    true_vertex[energy_index][1].pop_back();
	    true_vertex[energy_index][2].pop_back();
	    true_vertex[energy_index][3].pop_back();
	    reco_vertex[energy_index][0].pop_back();
	    reco_vertex[energy_index][1].pop_back();
	    reco_vertex[energy_index][2].pop_back();
	    reco_vertex[energy_index][3].pop_back();
	    error[energy_index][0].pop_back();
	    error[energy_index][1].pop_back();
	    //keep track of the number of events left out
	    events_excluded[energy_index]++;
	    events_per_file[x/5-1][y]--;
	  }
	  events_per_file[x/5-1][y]++;
	}
	myfile.close();	
      } else { 
	std::cout << "Unable to open file: " << filename << std::endl;
      }
    }
  }
  
  // --------------------------------------
  // Initialize zone sizes
  // ----------------------------------


  const int num_rzones = 16;
  const int num_zzones = 16;

  const double rstep = 21.875;  // = 350/num_rzones
  const double zstep = 50;     // = 800/num_zzones

  const int min_event_cutoff=150;
  
  std::cout << "num r zones: " << num_rzones << std::endl;
  std::cout << "num z zones: " << num_zzones << std::endl;
  std::cout << "min events in bin: expected amount" << std::endl;
  std::vector<double> rvals(num_rzones,0);
 for (int i=0;i<num_rzones;i++) {
   rvals.at(i)=rstep/2+(rstep*i);
 }
 std::vector<double> zvals(num_zzones,0);
 for (int i=0;i<num_zzones;i++) {
   zvals.at(i)=-400+(zstep/2)+(zstep*i);
 }


  // ---------------------------------------------
  // Print number of events per energy level
  // ---------------------------------------------
  for (int i=0; i<10; i++) {
    std::cout << i << " " << tb[i][0].size() << " " << reco_vertex[i][0].size() << std::endl;
  } 
 
 // ----------------------------------------------
 // Find the mean reco v error in each zone
 // ----------------------------------------------
 std::vector<std::vector<std::vector<double>>> init_mean_v_error(10, std::vector<std::vector<double>>(num_rzones, std::vector<double>(num_zzones,0)));
 // make 3D vector to keep track of the number of events in each z,r zone
 std::vector<std::vector<std::vector<double>>> count_init_events_in_zone_v_error(10, std::vector<std::vector<double>>(num_rzones, std::vector<double>(num_zzones,0)));
 // make 3D vector to keep track of which events  we add to which zone (for calculating standard deviation)
 // loop through every energy
 for (int energy=0; energy<10; energy++) {
   // loop through all the zones
   for (int r=0; r<357; r+=rstep) {
     int rindex = r/rstep;
     for (int z=-400; z<400; z+=zstep) {
       int zindex = (z+400)/zstep;
       // loop through every event and see which are in the current zone
       int length = reco_vertex[energy][3].size();
       for (int i=0; i<length; i++) {
         if (reco_vertex[energy][3].at(i) >= r && reco_vertex[energy][3].at(i) < r+rstep) {
           if (reco_vertex[energy][2].at(i) >= z && reco_vertex[energy][2].at(i) < z+zstep) {
             //if its in this zone, add it to that zone's mean r error and keep track of the number of events counted
             init_mean_v_error[energy][rindex][zindex] += sqrt(pow(reco_vertex[energy][0].at(i)-true_vertex[energy][0].at(i),2) + pow(reco_vertex[energy][1].at(i)-true_vertex[energy][1].at(i),2) + pow(reco_vertex[energy][2].at(i)-true_vertex[energy][2].at(i),2));
             count_init_events_in_zone_v_error[energy][rindex][zindex]++;
           }
         }
       } //end of event loop
       //before moving to the next zone, find the mean error of the current zone
       if (count_init_events_in_zone_v_error[energy][rindex][zindex] != 0) {
         init_mean_v_error[energy][rindex][zindex] = init_mean_v_error[energy][rindex][zindex]/count_init_events_in_zone_v_error[energy][rindex][zindex];
       }
     } //end of z zone loop
   } //end or r zone loop
 } //end of energy loop

 //-------------------------------
 //Find total number of events counted
 //--------------------------------
 std::cout << "checking... " << std::endl;
 std::vector<double> total_num_events(10,0);
 for (int e=0;e<10;e++) {
 	for (int r=0; r<num_rzones;r++) {
		for (int z=0;z<num_zzones;z++) {
			total_num_events[e]+=count_init_events_in_zone_v_error[e][r][z];
		}
	}
	std::cout << total_num_events[e] << std::endl;	
 }

 // ----------------------------------------------
 // Make 2D reco verto histo of error in each zone
 // ----------------------------------------------
 TH2D *init_verror5mev = new TH2D("init_verror5mev", "blank777", num_rzones, 0, 357, num_zzones, -400, 400);

 TH2D *init_verror50mev = new TH2D("init_verror50mev", "blank7777", num_rzones, 0, 357, num_zzones, -400, 400);

 for (int r=0; r<num_rzones; r++) {
   int cutoff0 = (2*r+1)*total_num_events[0]/3600;
   //make sure min cutoff is greater than min_event_cutoff as to not overcorrect
   if (cutoff0<min_event_cutoff) {
     cutoff0=min_event_cutoff;
   }
   int cutoff3 = (2*r+1)*total_num_events[3]/3600;
   if (cutoff3<min_event_cutoff) {
     cutoff3=min_event_cutoff;
   }
   for (int z=0; z<num_zzones; z++) {
     if (count_init_events_in_zone_v_error[0][r][z]>cutoff0) {
       init_verror5mev->Fill(rvals[r],zvals[z],init_mean_v_error[0][r][z]);
     } else {
       init_verror5mev->Fill(rvals[r],zvals[z],0);
     }
     if (count_init_events_in_zone_v_error[3][r][z]>cutoff3) {
       init_verror50mev->Fill(rvals[r],zvals[z],init_mean_v_error[3][r][z]);
     } else {
       init_verror50mev->Fill(rvals[r],zvals[z],0);
     }
   }
 }
 init_verror5mev->SetTitle("Initial Mean Error of Reconstructed Vertex at 5 MeV; Reconstructed Vertex R (mm); Reconstructed Vertex Z (mm); Mean Error of Reconstructed Vertex");
 init_verror5mev->Write();
 init_verror50mev->SetTitle("Initial Mean Error of Reconstructed Vertex at 20 MeV; Reconstructed Vertex R (mm); Reconstructed Vertex Z (mm); Mean Error of Reconstructed Vertex (mm)");
 init_verror50mev->Write();

  // ----------------------------------------------------
  // CALCULATE THE MEAN Z ERROR AS A FUNCTION OF R AND Z 
  // ----------------------------------------------------
  
  std::vector<std::vector<std::vector<double>>> mean_z_error(10, std::vector<std::vector<double>>(num_rzones, std::vector<double>(num_zzones,0)));
  // make 3D vector to keep track of the number of events in each z,r zone
  std::vector<std::vector<std::vector<double>>> count_events_in_zone_z_error(10, std::vector<std::vector<double>>(num_rzones, std::vector<double>(num_zzones,0))); 
  // make 3D vector to keep track of which events we add to which zone (for calculating standard deviation)
  std::vector<double> list_z_error[10][num_rzones][num_zzones]; 
  
  
  // loop through every energy
  for (int energy=0; energy<10; energy++) {
    // loop through all the zones
    for (int r=0; r<357; r+=rstep) {
      int rindex = r/rstep;
      for (int z=-400; z<400; z+=zstep) {
        int zindex = (z+400)/zstep;
	// loop through every event and see which are in the current zone
        int size = reco_vertex[energy][3].size();
	for (int i=0; i<size; i++) {
	  if (reco_vertex[energy][3].at(i) > r && reco_vertex[energy][3].at(i) < r+rstep) {
	    if (reco_vertex[energy][2].at(i) > z && reco_vertex[energy][2].at(i) < z+zstep) {
              //if its in this zone, add it to that zone's mean z error and keep track of the number of hits counted
              mean_z_error[energy][rindex][zindex] += error[energy][1].at(i);
	      count_events_in_zone_z_error[energy][rindex][zindex]++;
	      //also keep track of the value added for calculating standard deviation later
              list_z_error[energy][rindex][zindex].push_back( error[energy][1].at(i) );
	    } 
	  }
	} //end of event loop
	//before moving to the next zone, find the mean error of the current zone
	if (count_events_in_zone_z_error[energy][rindex][zindex] != 0) {
	  mean_z_error[energy][rindex][zindex] = mean_z_error[energy][rindex][zindex]/count_events_in_zone_z_error[energy][rindex][zindex];
        }
      } //end of z zone loop
    } //end or r zone loop
  } //end of energy loop

  // ----------------------------------------------------
  // CALCULATE THE STANDARD DEVIATION OF THE MEAN Z ERROR 
  // ----------------------------------------------------
  // now find the standard deviation of the mean z error at every energy
  std::vector<std::vector<std::vector<double>>> sd_z_error(10, std::vector<std::vector<double>>(num_rzones, std::vector<double>(num_zzones,0)));
  for (int energy=0; energy<10; energy++) {
    //loop through all the zones in our 7x8 grid
    for (int r=0; r<357; r+=rstep) {
      int rindex =r/rstep;
      for (int z=-400; z<400; z+=zstep) {
	int zindex = (z+400)/zstep;
        //loop through every event in that zone and sum the difference squared
	int size = count_events_in_zone_z_error[energy][rindex][zindex];
	for (int i=0; i<size; i++) {
          sd_z_error[energy][rindex][zindex] += pow( list_z_error[energy][rindex][zindex].at(i) - mean_z_error[energy][rindex][zindex], 2);
	}
	// divide by total number of entries in that zone and take the sqrt
	if (size !=0) {
	  sd_z_error[energy][rindex][zindex] = sd_z_error[energy][rindex][zindex]/size;
	  sd_z_error[energy][rindex][zindex] = sqrt( sd_z_error[energy][rindex][zindex] );
        }
      } //end of z zone loop
    } // end of r zone loop
  } // end of energy loop

 //--------------
 //Find total
 //--------------
  std::cout << "Checking again... " << std::endl;
  std::vector<double> tot_events (10,0);
  for (int e=0;e<10;e++) {
  	for (int r=0;r<num_rzones;r++) {
	     for (int z=0;z<num_zzones;z++) {
	     	  tot_events[e]+=count_events_in_zone_z_error[e][r][z];
	     }
	}
	std::cout << tot_events[e] << std::endl;
  }


 // -------------------------------------
 // Make 2D normal histo of zerror(recoz)
 // -------------------------------------- 
 TCanvas *c1 = new TCanvas("c1","correction graphs", 2000, 2000);
 c1->SetGrid();
 
 std::stringstream hzerror5mev_name;
 hzerror5mev_name << "hzerror5mev";
 TH2D hzerror5mev( hzerror5mev_name.str().c_str(), "Asymmetry Z Error; RecoZ; RecoZ - TrueZ", 800, -400, 400 , 800, -400, 400);

 int size = reco_vertex[0][2].size();
 for (int i=0; i<size; i++) {

   hzerror5mev.Fill(reco_vertex[0][2].at(i), error[0][1].at(i));
 }
 hzerror5mev.Write();
 
 std::stringstream hzerror50mev_name;
 hzerror50mev_name << "hzerror50mev";
 TH2D hzerror50mev( hzerror50mev_name.str().c_str(), "Asymmetry Z Error; RecoZ; RecoZ - TrueZ", 800, -400, 400 , 800, -400, 400);

 size = reco_vertex[9][2].size();
 for (int i=0; i<size; i++) {
   hzerror50mev.Fill(reco_vertex[9][2].at(i), error[9][1].at(i));
 }
 hzerror50mev.Write();

 // -----------------------------------------------------
 // Make 2D colz histo of the number of events in each zone
 // -----------------------------------------------------
 
 TH2D *hzr_numevents5mev = new TH2D("numevents5mev","blank1", num_rzones, 0, 357, num_zzones, -400, 400);

 TH2D *hzr_numevents50mev = new TH2D("numevents50mev","blank11", num_rzones, 0, 357, num_zzones, -400, 400);

 for (int r=0; r<num_rzones; r++) {
   for (int z=0; z<num_zzones; z++) {
     hzr_numevents5mev->Fill(rvals[r], zvals[z], count_events_in_zone_z_error[0][r][z]);
     hzr_numevents50mev->Fill(rvals[r], zvals[z], count_events_in_zone_z_error[3][r][z]);
   }
 }

 hzr_numevents5mev->SetTitle("Positions of Reconstructed Vertex at 5 MeV over 50,000 Events; Reconstructed Vertex R (mm); Reconstructed Vertex Z (mm); Number of Events");
 hzr_numevents5mev->Write();

 hzr_numevents50mev->SetTitle("Positions of Reconstructed Vertex at 20 MeV; Reconstructed R (mm); Reconstructed Z (mm); Number of Events");
 hzr_numevents50mev->Write(); 

 // -----------------------------------------
 // Make 2D colz histo of zerror(recor,recoz)
 // -----------------------------------------
 TH2D *hzr_zerror5mev = new TH2D("zerror5mev", "blank2", num_rzones, 0, 357, num_zzones, -400, 400);

 TH2D *hzr_zerror50mev = new TH2D("zerror50mev", "blank22", num_rzones, 0, 357, num_zzones, -400, 400);
 
 for (int r=0; r<num_rzones; r++) {
   for (int z=0; z<num_zzones; z++) {
     int cutoff0 = (2*r+1)*total_num_events[0]/3600;
     //make sure min cutoff is greater than min_event_cutoff as to not overcorrect
     if (cutoff0<min_event_cutoff) {
       cutoff0=min_event_cutoff;
     }
     int cutoff3 = (2*r+1)*total_num_events[3]/3600;
     if (cutoff3<min_event_cutoff) {
       cutoff3=min_event_cutoff;
     }
     if (count_events_in_zone_z_error[0][r][z]>cutoff0) { 
       hzr_zerror5mev->Fill(rvals[r],zvals[z],mean_z_error[0][r][z]);
     } else {
       hzr_zerror5mev->Fill(rvals[r],zvals[z],0);
     }
     if (count_events_in_zone_z_error[3][r][z]>cutoff3) {
       hzr_zerror50mev->Fill(rvals[r],zvals[z],mean_z_error[3][r][z]);
     } else {
       hzr_zerror50mev->Fill(rvals[r],zvals[z],0);
     }
   }
 }
 hzr_zerror5mev->SetTitle("Mean Error of Reconstructed Z Coordinate of Vertex at 5 MeV over 50,000 Events; Reconstructed Vertex R (mm); Reconstructed Vertex Z (mm); Mean Error of Reconstructed Vertex Z");
 hzr_zerror5mev->Write();
 hzr_zerror50mev->SetTitle("Mean Z Coordinate Error of Reconstructed Vertex at 20 MeV; Reconstructed R (mm); Reconstructed Z (mm); Mean Z Coordinate Error (mm)");
 hzr_zerror50mev->Write();
 
 // ------------------------------------------------------------------------------------
 // Correct results for all energies based on Z error ONLY IF THE BIN HAS ENOUGH RESULTS
 // ------------------------------------------------------------------------------------
 
 // loop through every energy
 for (int energy=0; energy<10; energy++) {
   // loop through all the zones
   for (int r=0; r<357; r+=rstep) {
     int rindex = r/rstep;
     for (int z=-400; z<400; z+=zstep) {
       int zindex = (z+400)/zstep;
       // loop through every event
       size = reco_vertex[energy][3].size();
       for (int i=0; i<size; i++) {
	 // see which are in current zone
         if (reco_vertex[energy][3].at(i) > r && reco_vertex[energy][3].at(i) < r+rstep) {
           if (reco_vertex[energy][2].at(i) > z && reco_vertex[energy][2].at(i) < z+zstep) {
             //Only correct it if there are enough events in the bin
	     int cutoff = ((2*rindex+1)*tot_events[energy])/3600;
	     if (cutoff<min_event_cutoff) {
  	       cutoff=min_event_cutoff;
	     }
	     if (count_events_in_zone_z_error[energy][rindex][zindex]>cutoff) {
               reco_vertex[energy][2].at(i)-=mean_z_error[energy][rindex][zindex];
 	     }
	   }
         }
       } //end of event loop 
     } //end of z zone loop
   } //end or r zone loop
 } //end of energy loop
 // ----------------------------------------------------
 // CALCULATE THE MEAN R ERROR AS A FUNCTION OF R AND Z 
 // ----------------------------------------------------
 std::vector<std::vector<std::vector<double>>> mean_r_error(10, std::vector<std::vector<double>>(num_rzones, std::vector<double>(num_zzones,0)));
 // make 3D vector to keep track of the number of events in each z,r zone
 std::vector<std::vector<std::vector<double>>> count_events_in_zone_r_error(10, std::vector<std::vector<double>>(num_rzones, std::vector<double>(num_zzones,0)));
 // make 3D vector to keep track of which events we add to which zone (for calculating standard deviation)
 // loop through every energy
 for (int energy=0; energy<10; energy++) {
   // loop through all the zones
   for (int r=0; r<357; r+=rstep) {
     int rindex = r/rstep;
     for (int z=-400; z<400; z+=zstep) {
       int zindex = (z+400)/zstep;
       // loop through every event and see which are in the current zone
       size = reco_vertex[energy][3].size();
       for (int i=0; i<size; i++) {
         if (reco_vertex[energy][3].at(i) > r && reco_vertex[energy][3].at(i) < r+rstep) {
           if (reco_vertex[energy][2].at(i) > z && reco_vertex[energy][2].at(i) < z+zstep) {
             //if its in this zone, add it to that zone's mean r error and keep track of the number of hits counted
             mean_r_error[energy][rindex][zindex] += error[energy][0].at(i);
             count_events_in_zone_r_error[energy][rindex][zindex]++;
           }
         }
       } //end of event loop 
       //before moving to the next zone, find the mean error of the current zone
       if (count_events_in_zone_r_error[energy][rindex][zindex] != 0) {
         mean_r_error[energy][rindex][zindex] = mean_r_error[energy][rindex][zindex]/count_events_in_zone_r_error[energy][rindex][zindex];
       }
     } //end of z zone loop
   } //end or r zone loop
 } //end of energy loop
 
 // -----------------------------------------------------
 // Make 2D colz histo of the number of events in each zone
 // -----------------------------------------------------
 TH2D *PZChzr_numevents5mev = new TH2D("PZCnumevents5mev","blank3", num_rzones, 0, 357, num_zzones, -400, 400);

 TH2D *PZChzr_numevents50mev = new TH2D("PZCnumevents50mev","blank33", num_rzones, 0, 357, num_zzones, -400, 400);

 for (int r=0; r<num_rzones; r++) {
   for (int z=0; z<num_zzones; z++) {
     PZChzr_numevents5mev->Fill(rvals[r], zvals[z], count_events_in_zone_r_error[0][r][z]);
     PZChzr_numevents50mev->Fill(rvals[r], zvals[z], count_events_in_zone_r_error[3][r][z]);
   }
 }
 PZChzr_numevents5mev->SetTitle("Post Z Correction Positions of Reconstructed Vertex at 5 MeV over 50,000 Events; Reconstructed Vertex R (mm); Reconstructed Vertex Z (mm); Number of Events");
 PZChzr_numevents5mev->Write();
 PZChzr_numevents50mev->SetTitle("Post Z Correction Positions of Reconstructed Vertex at 20 MeV over 50,000 Events; Reconstructed Vertex R (mm); Reconstructed Vertex Z (mm); Number of Events");
 PZChzr_numevents50mev->Write();

 // -----------------------------------------
 // Make 2D colz histo of rerror(recor,recoz)
 // -----------------------------------------
 TH2D *hzr_rerror5mev = new TH2D("rerror5mev", "blank4", num_rzones, 0, 357, num_zzones, -400, 400);
 
 TH2D *hzr_rerror50mev = new TH2D("rerror50mev", "blank44", num_rzones, 0, 357, num_zzones, -400, 400);
 
 for (int r=0; r<num_rzones; r++) {
   for (int z=0; z<num_zzones; z++) {
     int cutoff0 = (2*r+1)*total_num_events[0]/3600;
     //make sure min cutoff is greater than min_event_cutoff as to not overcorrect
     if (cutoff0<min_event_cutoff) {
       cutoff0=min_event_cutoff;
     }
     int cutoff3 = (2*r+1)*total_num_events[3]/3600;
     if (cutoff3<min_event_cutoff) {
       cutoff3=min_event_cutoff;
     }
     if (count_events_in_zone_r_error[0][r][z]>cutoff0) {
       hzr_rerror5mev->Fill(rvals[r],zvals[z],mean_r_error[0][r][z]);
     } else {
       hzr_rerror5mev->Fill(rvals[r],zvals[z],0);
     }
     if (count_events_in_zone_r_error[3][r][z]>cutoff3) {
       hzr_rerror50mev->Fill(rvals[r],zvals[z],mean_r_error[3][r][z]);
     } else {
       hzr_rerror50mev->Fill(rvals[r],zvals[z],0);
     }
   }
 }
 hzr_rerror5mev->SetTitle("Post Z Correction Mean Error of Reconstructed R Coordinate of Vertex at 50 MeV over 50,000 Events; Reconstructed Vertex R (mm); Reconstructed Vertex Z (mm); Mean Error of Reconstructed Vertex R");
 hzr_rerror5mev->Write();
 hzr_rerror50mev->SetTitle("Mean R Coordinate Error of Reconstructed Vertex at 20 MeV; Reconstructed R (mm); Reconstructed Z (mm); Mean R Coordinate Error (mm)");
 hzr_rerror50mev->Write();

 // -------------------------------------------------
 // Correct results for all energies based on R error 
 // -------------------------------------------------
 // loop through every energy
 for (int energy=0; energy<10; energy++) {
   // loop through all the zones 
   for (int r=0; r<357; r+=rstep) {
     int rindex = r/rstep;
     for (int z=-400; z<400; z+=zstep) {
       int zindex = (z+400)/zstep;
       // loop through every events
       size = reco_vertex[energy][3].size();
       for (int i=0; i<size; i++) {
	 // see which are in the current zone
         if (reco_vertex[energy][3].at(i) > r && reco_vertex[energy][3].at(i) < r+rstep) {
           if (reco_vertex[energy][2].at(i) > z && reco_vertex[energy][2].at(i) < z+zstep) {
             // only correct it if that zone saw enough events
	     int cutoff = ((2*rindex+1)*tot_events[energy])/3600;
	     if (cutoff<min_event_cutoff) {
  	       cutoff=min_event_cutoff;
	     }
	     if (count_events_in_zone_r_error[energy][rindex][zindex]>cutoff) {
               reco_vertex[energy][3].at(i)-=mean_r_error[energy][rindex][zindex];
             }
	   }
         }
       } //end of event loop 
     } //end of z zone loop
   } //end or r zone loop
 } //end of energy loop

  // -----------------------------------------------------
 // Make 2D colz histo of the number of events in each zone
 // -----------------------------------------------------
 TH2D *PRChzr_numevents5mev = new TH2D("PRCnumevents5mev","blank9", num_rzones, 0, 357, num_zzones, -400, 400);

 TH2D *PRChzr_numevents50mev = new TH2D("PRCnumevents50mev","blank99", num_rzones, 0, 357, num_zzones, -400, 400);

 for (int r=0; r<num_rzones; r++) {
   for (int z=0; z<num_zzones; z++) {
     PRChzr_numevents5mev->Fill(rvals[r], zvals[z], count_events_in_zone_r_error[0][r][z]);
     PRChzr_numevents50mev->Fill(rvals[r], zvals[z], count_events_in_zone_r_error[3][r][z]);
   }
 }
 PRChzr_numevents5mev->SetTitle("Post Z Correction Positions of Reconstructed Vertex at 5 MeV over 50,000 Events; Reconstructed Vertex R (mm); Reconstructed Vertex Z (mm); Number of Events");
 PRChzr_numevents5mev->Write();
 PRChzr_numevents50mev->SetTitle("Post Correction Positions of Reconstructed Vertex at 20 MeV; Reconstructed R (mm); Reconstructed Z (mm); Number of Events");
 PRChzr_numevents50mev->Write();


 // -------------------------------------------------------------------
 // Calculate new Reco X and Reco Y using the updated Reco R
 // -------------------------------------------------------------------
 //find theta of r guess for all events
 for (int e=0;e<10;e++) {
   size=reco_vertex[e][3].size();
   for (int i=0;i<size;i++) {
     double theta=atan2(reco_vertex[e][1].at(i),reco_vertex[e][0].at(i));
     reco_vertex[e][0].at(i)=reco_vertex[e][3].at(i)*cos(theta);
     reco_vertex[e][1].at(i)=reco_vertex[e][3].at(i)*sin(theta);
   }
 }

 // ----------------------------------------------
 // Find the mean reco v error in each zone
 // ----------------------------------------------
 std::vector<std::vector<std::vector<double>>> mean_v_error(10, std::vector<std::vector<double>>(num_rzones, std::vector<double>(num_zzones,0)));
 // make 3D vector to keep track of the number of events in each z,r zone
 std::vector<std::vector<std::vector<double>>> count_events_in_zone_v_error(10, std::vector<std::vector<double>>(num_rzones, std::vector<double>(num_zzones,0)));
 // make 3D vector to keep track of which events we add to which zone (for calculating standard deviation)
 // loop through every energy
 for (int energy=0; energy<10; energy++) {
   // loop through all the zones
   for (int r=0; r<357; r+=rstep) {
     int rindex = r/rstep;
     for (int z=-400; z<400; z+=zstep) {
       int zindex = (z+400)/zstep;
       // loop through every event and see which are in the current zone
       size = reco_vertex[energy][3].size();
       for (int i=0; i<size; i++) {
         if (reco_vertex[energy][3].at(i) > r && reco_vertex[energy][3].at(i) < r+rstep) {
           if (reco_vertex[energy][2].at(i) > z && reco_vertex[energy][2].at(i) < z+zstep) {
             //if its in this zone, add it to that zone's mean r error and keep track of the number of hits counted
             mean_v_error[energy][rindex][zindex] += sqrt(pow(reco_vertex[energy][0].at(i)-true_vertex[energy][0].at(i),2) + pow(reco_vertex[energy][1].at(i)-true_vertex[energy][1].at(i),2) + pow(reco_vertex[energy][2].at(i)-true_vertex[energy][2].at(i),2));
             count_events_in_zone_v_error[energy][rindex][zindex]++;
           }
         }
       } //end of events loop 
       //before moving to the next zone, find the mean error of the current zone
       if (count_events_in_zone_v_error[energy][rindex][zindex] != 0) {
         mean_v_error[energy][rindex][zindex] = mean_v_error[energy][rindex][zindex]/count_events_in_zone_v_error[energy][rindex][zindex];
       }
     } //end of z zone loop
   } //end or r zone loop
 } //end of energy loop

 // ----------------------------------------------
 // Make 2D reco verto histo of error in each zone
 // ----------------------------------------------
 TH2D *hzr_verror5mev = new TH2D("verror5mev", "blank7", num_rzones, 0, 357, num_zzones, -400, 400);

 TH2D *hzr_verror50mev = new TH2D("verror50mev", "blank77", num_rzones, 0, 357, num_zzones, -400, 400);

 for (int r=0; r<num_rzones; r++) {
   for (int z=0; z<num_zzones; z++) {
     int cutoff0 = (2*r+1)*total_num_events[0]/3600;
     //make sure min cutoff is greater than min_event_cutoff as to not overcorrect
     if (cutoff0<min_event_cutoff) {
       cutoff0=min_event_cutoff;
     }
     int cutoff3 = (2*r+1)*total_num_events[3]/3600;
     if (cutoff3<min_event_cutoff) {
       cutoff3=min_event_cutoff;
     }
     if (count_events_in_zone_v_error[0][r][z]>cutoff0) {
       hzr_verror5mev->Fill(rvals[r],zvals[z],mean_v_error[0][r][z]);
     } else {
       hzr_verror5mev->Fill(rvals[r],zvals[z],0);
     }
     if (count_events_in_zone_v_error[3][r][z]>cutoff3) {
       hzr_verror50mev->Fill(rvals[r],zvals[z],mean_v_error[3][r][z]);
     } else {
       hzr_verror50mev->Fill(rvals[r],zvals[z],0);
     }
   }
 }
 hzr_verror5mev->SetTitle("Post Z & R Correction Mean Error of Reconstructed Vertex at 50 MeV over 50,000 Events; Reconstructed Vertex R (mm); Reconstructed Vertex Z (mm); Mean Error of Reconstructed Vertex");
 hzr_verror5mev->Write();
 hzr_verror50mev->SetTitle("Post Correction Mean Error of Reconstructed Vertex at 20 MeV; Reconstructed R (mm); Reconstructed Z (mm); Mean Error of Reconstructed Vertex (mm)");
 hzr_verror50mev->Write();

 // -----------------------------------------------------
 // Make 2D colz histo of the number of events in each zone
 // -----------------------------------------------------
 TH2D *PZRChzr_numevents5mev = new TH2D("PZRCnumevents5mev","blank69", num_rzones, 0, 357, num_zzones, -400, 400);

 TH2D *PZRChzr_numevents50mev = new TH2D("PZRCnumevents50mev","blank699", num_rzones, 0, 357, num_zzones, -400, 400);

 for (int r=0; r<num_rzones; r++) {
   for (int z=0; z<num_zzones; z++) {
     PZRChzr_numevents5mev->Fill(rvals[r], zvals[z], count_events_in_zone_v_error[0][r][z]);
     PZRChzr_numevents50mev->Fill(rvals[r], zvals[z], count_events_in_zone_v_error[3][r][z]);
   }
 }
 PZRChzr_numevents5mev->SetTitle("Post Z & R Correction Positions of Reconstructed Vertex at 5 MeV over 50,000 Events; Reconstructed Vertex R (mm); Reconstructed Vertex Z (mm); Number of Events");
 PZRChzr_numevents5mev->Write();
 PZRChzr_numevents50mev->SetTitle("Post Correction Positions of Reconstructed Vertex at 20 MeV; Reconstructed R (mm); Reconstructed Z (mm); Number of Events");
 PZRChzr_numevents50mev->Write();

 // -----------------------------------------------------
 // Make 1D histo of the occurance of error at 5 and 50 mev
 // -----------------------------------------------------
 

 int ee=0;
 const int n =reco_vertex[ee][2].size();
 TH1D *herror5mev = new TH1D("error5mev","blank8", n, 0, 300);
 for (int i=0; i<n; i++) {
   herror5mev->Fill(sqrt(pow(reco_vertex[ee][0].at(i)-true_vertex[ee][0].at(i),2) + pow(reco_vertex[ee][1].at(i)-true_vertex[ee][1].at(i),2) + pow(reco_vertex[ee][2].at(i)-true_vertex[ee][2].at(i),2)));
 }

 
 
 ee=3;
 const int m = reco_vertex[ee][2].size();
 TH1D *herror50mev = new TH1D("error50mev","blank88", m, 0, 300);
 for (int i=0; i<m; i++) {
   herror50mev->Fill(sqrt(pow(reco_vertex[ee][0].at(i)-true_vertex[ee][0].at(i),2) + pow(reco_vertex[ee][1].at(i)-true_vertex[ee][1].at(i),2) + pow(reco_vertex[ee][2].at(i)-true_vertex[ee][2].at(i),2)));

 }

 herror5mev->SetTitle("Error of Reconstructed Vertex at 5 MeV");
 herror5mev->Write();
 herror50mev->SetTitle("Error of Reconstructed Vertex at 20 MeV");
 herror50mev->Write();






 // -----------------------------------------
 // Find the mean vertex error at each energy 
 // -----------------------------------------
 std::vector<double> c_vertex(10,0);

 std::vector<double> mean_vertex_error(10,0);
 //make a list of values used to make sd
 std::vector<double> list_vertex_error[10];
 // loop through every energy
 for (int e=0;e<10;e++) {
   // loop through every event at that energy
   size=reco_vertex[e][2].size();
   for (int i=0;i<size;i++) {
     // find the 3d vertex error and add it to the mean
     list_vertex_error[e].push_back(sqrt(pow(reco_vertex[e][0].at(i)-true_vertex[e][0].at(i),2) + pow(reco_vertex[e][1].at(i)-true_vertex[e][1].at(i),2) + pow(reco_vertex[e][2].at(i)-true_vertex[e][2].at(i),2)));
     mean_vertex_error[e]+=list_vertex_error[e].back();
   } //end of event loop 
   mean_vertex_error[e]=mean_vertex_error[e]/size;
   //find error which 68% of results are better than
   std::sort(list_vertex_error[e].begin(), list_vertex_error[e].end());
   c_vertex[e]=list_vertex_error[e].at(0.68*size);
 } //end of energy loop
 





  // ----------------------------------------------------
  // CALCULATE THE STANDARD DEVIATION OF THE MEAN VERTEX ERROR 
  // ----------------------------------------------------
  // now find the standard deviation of the mean vertex error at every energy
  std::vector<double> sd_vertex_error(10,0);
  for (int e=0;e<10;e++) {
    //loop through every event in that zone and sum the difference squared
    size=reco_vertex[e][2].size();
    for (int i=0;i<size;i++) {
      sd_vertex_error[e]+=pow(list_vertex_error[e].at(i)-mean_vertex_error[e],2);
    }
    // divide by total number of entries for that energy and take the sqrt
    sd_vertex_error[e]=sd_vertex_error[e]/size;
    sd_vertex_error[e]=sqrt(sd_vertex_error[e]);
  } // end of energy loop
  // -------------------------------------------------
  // Plot the mean vertex error as a function o energy
  // -------------------------------------------------
  TCanvas *c2 = new TCanvas("c2","reco vertex error", 5000, 5000);
  c2->cd(); 
  c2->SetGrid();

  //make energies to plot
  std::vector<double> energies = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50};
  std::vector<double> blank(10,0);
  
  for (int i=0;i<10;i++) {
    std::cout << energies.at(i) << " " << c_vertex.at(i) << " " << mean_vertex_error.at(i) << " " << blank.at(i) << " " << sd_vertex_error.at(i) << std::endl;
  }

  TGraphErrors *vertex_error = new TGraphErrors(10,energies.data(),mean_vertex_error.data(),blank.data(),sd_vertex_error.data());
  vertex_error->SetTitle("Mean Spatial Resolution; Energy of e- Produced (MeV); Mean Error of Reconstructed Vertex (mm)");
  vertex_error->SetMarkerColor(4);
  vertex_error->SetMarkerStyle(21);
  vertex_error->SetName("vertex_error");
  vertex_error->Write();

  TGraph *cvertex = new TGraph(10,energies.data(),c_vertex.data());
  cvertex->SetTitle("Absolute Spatial Resolution; Energy of e- Produced (MeV); Absolute Error of Reconstructed Vertex (mm)");
  cvertex->SetMarkerColor(4);
  cvertex->SetMarkerStyle(21);
  cvertex->SetName("cvertex");
  cvertex->Write();
  
  std::cout << " " << std::endl;
  std::cout << "events excluded (reco z<200 & z>-200): " << std::endl;
  for (int i=0;i<10;i++) {
    std::cout  << events_excluded.at(i) << std::endl;
  } 

  fout->Close();

  // ------------------------------------------------
  // Create text files to write corrected reconstructed vertex to
  // ------------------------------------------------
   
  std::cout << "final check... " << std::endl;
  for (int e=0;e<10;e++) {
     std::cout << reco_vertex[e][0].size() << std::endl;
  }



  for (int x=0;x<10;x++) {
      std::string filename = "/cluster/tufts/wongjiradlabnu/bharri06/data/corr_data/data_";
      filename += std::to_string(5*x+5);
      filename += "mev.txt";
      //std::cout << "Filename: " << filename << std::endl;
      std::ofstream outputfile;
      outputfile.open(filename, std::ios::trunc);
      //cout mean z error to a file
      outputfile << "Mean Z Error " << 5*x+5 << " MeV" << std::endl;
      for (int z=0;z<num_zzones;z++) {
	for (int r=0;r<num_rzones;r++) {
          //only if that zone saw enough events
          int cutoff = ((2*r+1)*tot_events[x])/3600;
          if (cutoff<min_event_cutoff) {
               cutoff=min_event_cutoff;
          }
	  if (count_events_in_zone_z_error[x][r][z]>cutoff) {
	    outputfile << mean_z_error[x][r][z] << " ";
	  } else {
            outputfile << 0 << " ";
	  } 
	}
	outputfile << std::endl;
      }
      outputfile << std::endl;
      //cout mean r error to a file
      outputfile << "Mean R Error " << 5*x+5 << " MeV" << std::endl;
      for (int z=0;z<num_zzones;z++) { 
        for (int r=0;r<num_rzones;r++) {
          //only if that zone saw enough events
          int cutoff = ((2*r+1)*tot_events[x])/3600;
          if (cutoff<min_event_cutoff) {
               cutoff=min_event_cutoff;
          }
          if (count_events_in_zone_r_error[x][r][z]>cutoff) {
            outputfile << mean_r_error[x][r][z] << " ";
          } else {
            outputfile << 0 << " ";
          }
        }
        outputfile << std::endl;
      }
      outputfile.close();
  }

  return 0;
}
