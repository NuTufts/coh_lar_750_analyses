// c++ headers
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
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
#include "/cluster/tufts/wongjiradlabnu/bharri06/cenns10geant4/analysis/include/CENNSEvent.hh"
#include "/cluster/tufts/wongjiradlabnu/bharri06/cenns10geant4/analysis/include/CENNSChannel.hh"

double* rotate_around_x(double x, double y, double z, double rotation);
double* rotate_around_y(double x, double y, double z, double rotation);
double* rotate_around_z(double x, double y, double z, double rotation);
double angle_between(double x1, double y1, double z1, double x2, double y2, double z2);

int main( int nargs, char** argv ) {
  // ---------------------------------------------------------------
  // Check to make sure the correct number of arguments are specified
  // ---------------------------------------------------------------
  std::cout << "Reconstructing the primary vertex..." << std::endl;
  if ( nargs!=3 ) {
    std::cout << "usage: reco_vertex [input cenns root file] [output root file]" << std::endl;
    return 0;
  }
  std::string input_cenns_geant_file = argv[1];
  std::string output_root = argv[2];
  
  TApplication app( "app", &nargs, argv );
  // -------------------------------------------------------
  // Open the input file and get branch data from each event
  // -------------------------------------------------------
  TFile* finput = new TFile( input_cenns_geant_file.c_str(), "OPEN" );
  //std::cout << "Input file opened." << std::endl;
  TTree* cenns = (TTree*)finput->Get("CENNS");
  //std::cout << "Finput pointed to Get(CENNS)" << std::endl;
  CENNSEvent* eventdata = nullptr;
  //std::cout << "Nullptr set to eventdata" << std::endl;
  TBranch* b_event = 0;
  //std::cout << "Nullptr set to branch event" << std::endl;
  cenns->SetBranchAddress("Event",&eventdata,&b_event);
  //std::cout << "Branch address set." << std::endl;
  int nentries = cenns->GetEntries();
  //std::cout << "Number of entries in file: " << nentries << std::endl;
  // ----------------------------------------
  // Get energy_batch# from input file name
  // ----------------------------------------
  std::cout << "Input file: " << input_cenns_geant_file << " (length = " << input_cenns_geant_file.size() << " )" << std::endl; 
  input_cenns_geant_file.erase(0, 57) ;
  //std::cout << "Input file: " << input_cenns_geant_file << " (length = " << input_cenns_geant_file.size() << " )" << std::endl;

  if (input_cenns_geant_file.size() == 15) {
  	input_cenns_geant_file.erase(1,7);
	input_cenns_geant_file.erase(3,5);
  } 
  if (input_cenns_geant_file.size() == 16) {
	if (input_cenns_geant_file.at(8) != '_') {
		input_cenns_geant_file.erase(2,7);
		input_cenns_geant_file.erase(4,5);
	} else {
		input_cenns_geant_file.erase(1,7);
		input_cenns_geant_file.erase(4,5);
	}  
  }
  if (input_cenns_geant_file.size() == 17) {
  	input_cenns_geant_file.erase(2,7);
	input_cenns_geant_file.erase(5,5);
  } 
  std::cout << "File: " << input_cenns_geant_file << std::endl;
  
  // ------------------------
  // CREATE OUTPUT ROOT FILE
  // ------------------------
  TFile* fout = new TFile( output_root.c_str(), "NEW" );
  
  // ---------------------------------------------------------------------------------------
  // Open a text file containing the postion of each PMT, and store the values in a 2D array
  // ---------------------------------------------------------------------------------------
  double pmt_center[114][3] = {0}; //top or bottom, square in 9x9 array, xyz coordinate
  std::ifstream pmtpositions ("/cluster/tufts/wongjiradlab/bharri06/cenns10geant4/gdml/lar1ton_pmt_and_sipm/pythonscripts/pmt_hole_positions.xml");
  if (pmtpositions.is_open()) {
     std::string temp_a, temp_b, x, y, z, temp_c, temp_d;
     int line_count = 0;
     while (pmtpositions >> temp_a >> temp_b >> x >> y >> z >> temp_c >> temp_d) {
         x.erase(x.begin(), x.begin()+3);
         x.pop_back();
         y.erase(y.begin(), y.begin()+3);
         y.pop_back();
         pmt_center[line_count][0] = std::stod(x);
         pmt_center[line_count][1] = std::stod(y);
         pmt_center[line_count][2] = 457.2;
         pmt_center[line_count+57][0] = std::stod(x);
         pmt_center[line_count+57][1] = std::stod(y);
         pmt_center[line_count+57][2] = -457.2;
         line_count++;
     }
  } else {
     std::cout << "not open" << std::endl;
  }
  pmtpositions.close();
  // ------------------------------------------------
  // Create text file to write reconstructed vertex to
  // ------------------------------------------------
  std::string filename = "/cluster/tufts/wongjiradlabnu/bharri06/data/reco_data/final2/data_";
  filename += input_cenns_geant_file;
  filename += ".txt";
  std::cout << "Filename: " << filename << std::endl;
  std::ofstream outputfile;
  outputfile.open(filename, std::ios::trunc);
  // ------------------------------------------------------------------
  // Erase batch# from input file, leaving just the energy of the event
  // ------------------------------------------------------------------
  if (input_cenns_geant_file.size() == 3) {
  	input_cenns_geant_file.erase(input_cenns_geant_file.begin()+1, input_cenns_geant_file.end());
  } else if (input_cenns_geant_file.size() == 4) {
	if (input_cenns_geant_file.at(1) != '_') {
		input_cenns_geant_file.erase(input_cenns_geant_file.begin()+2, input_cenns_geant_file.end());
	} else {
		input_cenns_geant_file.erase(input_cenns_geant_file.begin()+1, input_cenns_geant_file.end());
	}
  } else if (input_cenns_geant_file.size() == 5) {
  	input_cenns_geant_file.erase(input_cenns_geant_file.begin()+2, input_cenns_geant_file.end());
  }
  int energy = std::stoi(input_cenns_geant_file);
  std::cout << "Energy: " << input_cenns_geant_file << std::endl;

  // -------------------------------------
  // Loop through every event in input file
  // -------------------------------------
  for ( int ientry=0; ientry < nentries; ientry++ ) {
    std::cout << "------------------------------------------------------------------------------" << std::endl;
    std::cout << "[ENTRY " << ientry << "] " << std::endl;
    ULong_t bytes = cenns->GetEntry(ientry);
    if ( bytes==0 ) {
      // empty event
      continue;
    }
    // -----------------------------------
    // Initialize event specific variables
    // -----------------------------------
    std::vector<double> ceren_hit_x;   //x location of each cherenkov hit (push back value)
    std::vector<double> ceren_hit_y;   //y location of each cherenkov hit (push back value)
    std::vector<double> ceren_hit_z;   //z location of each cherenkov hit (push back value)
    
    std::vector<double> hit[3];         //[x,y,z] location of each scintillation hit (push back value)
    double scint_hit_count[114] = {0};  //number of scint hits on each PMT (increase value each time - start at zero)
    double tb_all[2] = {0};                 //[tob, bottom] number of scint hits on each endcap (increase value each time - start at zero)

    std::vector<double> channel_tracker;

    std::stringstream hxzrot_true_name;
    hxzrot_true_name << "hxzrot_true_entry" << ientry;
    TH2D hxzrot_true( hxzrot_true_name.str().c_str(), "score; x rot; z rot", 181, 0, 180, 361, 0 , 360);
    std::stringstream hxzrot_reco_name;
    hxzrot_reco_name << "hxzrot_reco_entry" << ientry;
    TH2D hxzrot_reco( hxzrot_reco_name.str().c_str(), "score; x rot; z rot", 181, 0, 180, 361, 0 , 360);


    // -----------------------------------
    // Loop through every channel in event
    // ----------------------------------
    for ( auto const& chdata : eventdata->fvChannels ) {
      int channelid = chdata.GetChannelNum();
      // std::cout << "  Channel[" << channelid << "] numhits=" << numhits << std::endl;
      // ---------------------------------
      // Loop through every hit in channel
      // ---------------------------------
      for ( auto const& photonhit : chdata ) {
        // Required variables for QE of photosensor
        int is_ceren = photonhit.ischerenkov;
        // Get hit positions
        double HitPosX = photonhit.pos[0];
        double HitPosY = photonhit.pos[1];
        double HitPosZ = photonhit.pos[2];
        // Only use scintillation hits to reconstruct the vertex
        if (!is_ceren) {
          // store hit position
	  hit[0].push_back( HitPosX );
          hit[1].push_back( HitPosY );
          hit[2].push_back( HitPosZ );
	  // keep track of which PMT was hit
	  scint_hit_count[channelid]++;
	  // keep track of if hit was on the top of bottom
	  if ( HitPosZ > 0 ) {
            tb_all[0]++;
          } else {
            tb_all[1]++;
          }
	} else {
          ceren_hit_x.push_back( HitPosX );
          ceren_hit_y.push_back( HitPosY );
          ceren_hit_z.push_back( HitPosZ );
	  //keep track of which channel was hit
          channel_tracker.push_back(channelid);
	}
      }// end of hit loop
    }//end of channel loop
    // ----------------------------------------------------------------------------------------------
    // Write the energy, number of scint hits on the top and bottom, and true vertex to the outputfile
    // ----------------------------------------------------------------------------------------------
    double xtruth = eventdata->fGenX;
    double ytruth = eventdata->fGenY;
    double ztruth = eventdata->fGenZ;
    outputfile << energy << " " << tb_all[0]<< " " << tb_all[1] << " " << xtruth << " " << ytruth << " "  << ztruth << " ";
    // ------------------------------------------------------------------------------------------------------------------------------
    // Reconstruct the z coordinate of the vertex using the asymmetry of number of hits on the top and bottom endcaps of the detector
    // ------------------------------------------------------------------------------------------------------------------------------
    double reco_vertex[3] = {0, 0, 0};
    double asymmetry = (tb_all[0]-tb_all[1])/(tb_all[0]+tb_all[1]);
    //4th degree polynomial fit (not the best idea because its obviously linear)
    //reco_vertex[2] = 1.56991 + 685.298*asym + 2.04786*pow(asym,2) - 309.799*pow(asym,3) + 13.8992*pow(asym,4);
    //linear fit of data
    //make sure to not guess above 400 or below -400
    if (asymmetry > 0.65) {
      asymmetry=0.65;
    }
    if (asymmetry < -0.65) {
      asymmetry=-0.65;
    }
    reco_vertex[2] = 615.384615385*asymmetry;
    // ---------------------------------------------------------------------------------------------------------------------------------
    // Make initial guess of x and y coordinate of the vertex using the top four most hit PMT on the endcap which recieved the most hits
    // ---------------------------------------------------------------------------------------------------------------------------------
    const int pmt_kept = 4;
    //Find ^ most hit PMT on top and bottom endcap of detector
    int max[2][pmt_kept] = {{0, 0, 0, 0}, {0, 0, 0, 0}};
    int index[2][pmt_kept] = {{0, 0, 0, 0}, {0, 0, 0, 0}};
    for (int square = 0; square < 114; square++) {
      int tb;
      if (square < 57) {
        tb = 0;
      } else {
        tb = 1;
      }
      int size = scint_hit_count[square];
      bool been_placed = false;
      for (int pmt = 0; pmt < pmt_kept ; pmt++) {
        if (size > max[tb][pmt] && !been_placed) {
          for (int j = pmt_kept-1; j > pmt; j--) {
            max[tb][j] = max[tb][j-1];
            index[tb][j]= index[tb][j-1];
          }
          max[tb][pmt] = size;
          index[tb][pmt] = square;
          been_placed = true;
        }
      }
    }
    // Use the most hit PMT on the endcap which recieved the most hits to reconstruct x & y coordinates of vertex
    // Was written this way to test if should fit all data, or only data on the endcap which recieved the most hits
    int total_hits_tb[2] = {0, 0};
    for (int tb = 0; tb < 2; tb++) {
      for (int pmt = 0; pmt < pmt_kept; pmt++) {
        total_hits_tb[tb] = total_hits_tb[tb] + max[tb][pmt];
      }
    }
    int tb;
    if (total_hits_tb[0] > total_hits_tb[1]) {
      tb = 0;
    } else {
      tb = 1;
    }
    for (int xy = 0; xy < 2; xy++) {
      for (int pmt = 0; pmt < pmt_kept; pmt++) {
        for (int j = 0; j < max[tb][pmt]; j++) {
          reco_vertex[xy] = reco_vertex[xy] + pmt_center[index[tb][pmt]][xy];
        }
      }
      reco_vertex[xy] = reco_vertex[xy]/total_hits_tb[tb];
    }
    // ---------------------------------------------------------------------------------------------
    // Use that intial guess of the x and y coordinates to reconstruct again using gradient descent
    // ---------------------------------------------------------------------------------------------
    // Use the endcap with the most hits
    int start = 0;
    int end = 114;
    if (tb_all[0] > tb_all[1]) {
      end = 57;
    } else {
      start = 57;
    }
    // Initialize standard deviation variables
    double sigma_temp_x = 0;
    double sigma_temp_y = 0;
    double sigma_temp = 0;
    int num_added = 0;
    // Calculate standard deviation of initial x and y coordinate guesses
    for (int pmt=start; pmt<end; pmt++) {
      int size = scint_hit_count[pmt];
      for (int i=0; i<size; i++) {
        sigma_temp_x += pow(pmt_center[pmt][0]/357-reco_vertex[0]/357,2);  //  +pow(pmt_center[pmt][1]/357-true_vertex[1]/357,2);
        sigma_temp_y += pow(pmt_center[pmt][1]/357-reco_vertex[1]/357,2);
	num_added++;
      }
    }
    // Since the distribution is expected to be spherical, we average our x and y standard deviations to use as an overarching sigma in the loss function to create spherical results
    sigma_temp_x = sqrt(sigma_temp_x/num_added);
    sigma_temp_y = sqrt(sigma_temp_y/num_added);
    sigma_temp = (sigma_temp_x+sigma_temp_y)/2; 
    // Initialize the loss function variables
    double error = 0;
    double error_temp = 0;
    double sigma = 0;
    double a = 0;
    double gradient_sigma(0);
    double gradient_a(0);
    std::vector<double> hits_observed(57, 0);
    std::vector<double> mu(2,0);
    std::vector<double> hits_predicted(57, 0);
    std::vector<double> gradient_mu(2,0);
    //keep track of most hit pmt to standardize the amplitude, a
    double most_hit_pmt = 0;
    for (int pmtid=0; pmtid<57; pmtid++) {
      //standardize the number of hits observed
      if (tb_all[0] > tb_all[1]) {
        hits_observed[pmtid] = scint_hit_count[pmtid];
      } else {
	hits_observed[pmtid] = scint_hit_count[pmtid+57];
      }
      if (hits_observed[pmtid] > most_hit_pmt) {
	most_hit_pmt = hits_observed[pmtid];
      }
    }
    for (int pmtid=0; pmtid<57; pmtid++) {
      hits_observed[pmtid] = hits_observed[pmtid]/most_hit_pmt;
    }
    // make initial guess of loss function variables
    mu[0] = reco_vertex[0]/357;
    mu[1] = reco_vertex[1]/357;
    sigma = sigma_temp;
    a = 1;
    // Use gradient descent to minimize loss function of fitted 2D gausian-esq function
    // Find gradients of all the variables
    for (int pmtid=0; pmtid<57; pmtid++) {
      double base_exp_term = exp(-pow(pmt_center[pmtid][0]/357-mu[0],2)/(2*pow(sigma,2)))*exp(-pow(pmt_center[pmtid][1]/357-mu[1],2)/(2*pow(sigma,2)));
      hits_predicted[pmtid] = a*base_exp_term;
      error += pow(hits_observed[pmtid]-hits_predicted[pmtid],2);
      double base_grad_term = -2*(hits_observed[pmtid]-hits_predicted[pmtid]);
      gradient_a += base_grad_term*base_exp_term;
      gradient_sigma += base_grad_term*(pow(pmt_center[pmtid][0]/357-mu[0],2)+pow(pmt_center[pmtid][1]/357-mu[1],2))*hits_predicted[pmtid]/pow(sigma,3);
      gradient_mu[0] += base_grad_term*-hits_predicted[pmtid]*(pmt_center[pmtid][0]/357-mu[0])/pow(sigma,2);
      gradient_mu[1] += base_grad_term*-hits_predicted[pmtid]*(pmt_center[pmtid][1]/357-mu[1])/pow(sigma,2);
    }
    error_temp = error;
    // Set the learning rates for all the variables
    int iteration = 1;
    double lr_a = abs((0.001*a)/gradient_a);
    double lr_sigma = abs((0.001*sigma)/gradient_sigma);
    double lr_mu_x = abs((0.001*mu[0])/gradient_mu[0]);
    double lr_mu_y = abs((0.001*mu[1])/gradient_mu[1]);
    double lratemp = lr_a;
    double lrsigmatemp = lr_sigma;
    double lrmuxtemp = lr_mu_x;
    double lrmuytemp = lr_mu_y;
    // cout variables and their gradients before starting for debugging
    //std::cout << "Initial Values:" << std::endl;
    //std::cout << "               A: " << a << "   GA: " << gradient_a << std::endl;	
    //std::cout << "               S: " << sigma << "   GS: " << gradient_sigma << std::endl;
    //std::cout << "Top Hits: " << tb_all[0] << "  Bot Hits: " << tb_all[1] << std::endl;
    //std::cout << "TrY: " << eventdata->fGenY/357 << "  muY: " << mu[1] << "  GmuY: " << gradient_mu[1] <<  std::endl;
    //std::cout << "Error: " << error << std::endl;
    //std::cout << "lr_a: " << lr_a << "  lr_sig: " << lr_sigma << "  lr mu x: " << lr_mu_x << "  lr mu y: " << lr_mu_y << std::endl;
    std::vector<double> mu_tracker[2];
    mu_tracker[0].push_back(mu[0]);
    mu_tracker[1].push_back(mu[1]);
    // While gradient magnitude is less than desired value update all values if they're on
    bool on = true;
    double prev_change_a = 0;
    double prev_change_sigma = 0;
    std::vector<double> prev_change_mu(2,0);
    double momentum_a = 0;
    double momentum_sigma = 0;
    double momentum_mu = 0;
    double a_reset = 1;
    double sigma_reset = 1;
    std::vector<double> mu_reset(2,1);
    std::vector<double> best_coordinates = {99, 0, 0, 0, 0};
    // Do ten thousand iterations and pick the one with the lowest error 
    while (on && iteration < 10000 ) { 
    //while ( sqrt(pow(gradient_mu[time_cut][0],2)+pow(gradient_mu[time_cut][1],2)) > 1.0e-4 ) {
    //while ( sqrt(pow(gradient_A[time_cut],2)+pow(gradient_sigma[time_cut],2)+pow(gradient_mu[time_cut][0],2)+pow(gradient_mu[time_cut][1],2)) > 1.0e-4 ) {
      //double lr = 1.0e-3/(iteration+1); //(iteration+1)/500.0;
      // find which direction is best
      a -= lr_a*gradient_a+prev_change_a*momentum_a;
      sigma -= lr_sigma*gradient_sigma+prev_change_sigma*momentum_sigma;
      mu[0] -= lr_mu_x*gradient_mu[0]+prev_change_mu[0]*momentum_mu;
      mu[1] -= lr_mu_y*gradient_mu[1]+prev_change_mu[1]*momentum_mu;
      mu_tracker[0].push_back(mu[0]);
      mu_tracker[1].push_back(mu[1]);
      prev_change_a = lr_a*gradient_a+prev_change_a*momentum_a;
      prev_change_sigma = lr_sigma*gradient_sigma+prev_change_sigma*momentum_sigma;
      prev_change_mu[0] = lr_mu_x*gradient_mu[0]+prev_change_mu[0]*momentum_mu;
      prev_change_mu[1] = lr_mu_y*gradient_mu[1]+prev_change_mu[1]*momentum_mu;
      // make temp variables to hold previous gradients for new lr
      double temp_gradient_a = gradient_a;
      double temp_gradient_sigma = gradient_sigma;
      std::vector<double> temp_gradient_mu = {gradient_mu[0], gradient_mu[1]};
      // reset gradients
      error = 0;
      gradient_sigma = 0;
      gradient_mu[0] = 0;
      gradient_mu[1] = 0;
      // calculate new gradients
      for (int pmtid=0; pmtid<57; pmtid++) {
        double base_exp_term = exp(-pow(pmt_center[pmtid][0]/357-mu[0],2)/(2*pow(sigma,2)))*exp(-pow(pmt_center[pmtid][1]/357-mu[1],2)/(2*pow(sigma,2)));
        hits_predicted[pmtid] = a*base_exp_term;
        error += pow(hits_observed[pmtid]-hits_predicted[pmtid],2);
        double base_grad_term = -2*(hits_observed[pmtid]-hits_predicted[pmtid]);
        gradient_a += base_grad_term*base_exp_term;
        gradient_sigma += base_grad_term*(pow(pmt_center[pmtid][0]/357-mu[0],2)+pow(pmt_center[pmtid][1]/357-mu[1],2))*hits_predicted[pmtid]/pow(sigma,3);
        gradient_mu[0] += base_grad_term*hits_predicted[pmtid]*(pmt_center[pmtid][0]/357-mu[0])/pow(sigma,2);
        gradient_mu[1] += base_grad_term*hits_predicted[pmtid]*(pmt_center[pmtid][1]/357-mu[1])/pow(sigma,2);
      }
      // ---------------------------------------------------------------------------
      // Store the coordinates with the lowest error
      // Make sure the r coordinate is inside the detector!!!!!
      // ---------------------------------------------------------------------------
      if ( (error < best_coordinates[0]) && (357 > sqrt(pow(mu[0]*357,2)+pow(mu[1]*357,2))) ) {
        best_coordinates[0] = error;
	best_coordinates[1] = a;
	best_coordinates[2] = sigma;
	best_coordinates[3] = mu[0];
	best_coordinates[4] = mu[1];
      }
      // check to see if sign of gradient switched so we can decrease the learning rate
      if ((abs(gradient_a) < 1.0e-2) and gradient_a/temp_gradient_a < 0) {
	//std::cout << "------------------------a reset triggered----------------------" << std::endl;
	a_reset = a_reset*10;
      }
      if ((abs(gradient_sigma) < 1.0e-2) and gradient_sigma/temp_gradient_sigma < 0) {
	//std::cout << "------------------------sigma reset triggered-----------------" << std::endl;
	sigma_reset = sigma_reset*10;;
      }
      if ((abs(gradient_mu[0]) < 1.0e-2) and gradient_mu[0]/temp_gradient_mu[0] < 0) { 
	//std::cout << "------------------------mux reset triggered-----------------" << std::endl;
	mu_reset[0] = mu_reset[0]*10;
      }
      if ((abs(gradient_mu[1]) < 1.0e-2) and gradient_mu[1]/temp_gradient_mu[1] < 0) {
	//std::cout << "------------------------muy reset triggered-----------------" << std::endl;
	mu_reset[1] = mu_reset[1]*10;
      }
      if (error_temp*1.5 < error) {
	on = false;
      } else {
	error_temp = error;
      }
      // calculate new learning rates
      double i = iteration;
      double num1 = i/240;
      double num2 = i/80;
      double num3 = i/70;
      lr_a = lratemp*pow(2.71828,num3);
      lr_a = lr_a/a_reset;
      lr_sigma = lrsigmatemp*pow(2.71828,num2);
      lr_sigma = lr_sigma/sigma_reset;
      lr_mu_x = lrmuxtemp*pow(2.71828,num1);
      lr_mu_x = lr_mu_x/mu_reset[0];
      lr_mu_y = lrmuytemp*pow(2.71828,num1);
      lr_mu_y = lr_mu_y/mu_reset[1];
      // print variables after each iteration for debugging
      //std::cout << "exiting loop" << std::endl;
      //std::cout << "Values after " << iteration << " iterations" << std::endl;
      //std::cout << "Initial Values:" << std::endl;
      //std::cout << "               A: " << a << "   GA: " << gradient_a << std::endl;
      //std::cout << "               S: " << sigma << "   GS: " << gradient_sigma << std::endl;
      //std::cout << "TrX: " << eventdata->fGenX/357 << "  muX: " << mu[0] << "  GmuX: " << gradient_mu[0] << std::endl;
      //std::cout << "TrY: " << eventdata->fGenY/357 << "  muY: " << mu[1] << "  GmuY: " << gradient_mu[1] <<  std::endl;
      //std::cout << "Z: " << eventdata->fGenZ << "  Asym: " << asym << "  Top Hits: " << tb_all[0] << "  Bot Hits: " << tb_all[1] << std::endl;
      //std::cout << "Error: " << error << std::endl;
      //std::cout << "lr_a: " << lr_a << "  lr_sig: " << lr_sigma << "  lr mu x: " << lr_mu_x << "  lr mu y: " << lr_mu_y << std::endl;
      //std::cout << std::endl;
      iteration++;
    }
    // print variables after the final iteration
    //std::cout << "Values after " << iteration << " iterations" << std::endl;
    //std::cout << "               A: " << a << "   GA: " << gradient_a << std::endl;
    //std::cout << "               S: " << sigma << "   GS: " << gradient_sigma << std::endl;
    //std::cout << "TrX: " << eventdata->fGenX/357 << "  muX: " << mu[0] << "  GmuX: " << gradient_mu[0] << std::endl;
    
    // ---------------------------------------------------------------------------
    // Convert normalized mu_x and mu_y to x and y coordinates inside the detector
    // ---------------------------------------------------------------------------
    if (sqrt(pow(best_coordinates[3]*350,2)+pow(best_coordinates[4]*350,2)) > 350) {
      double theta = atan2(best_coordinates[4],best_coordinates[3]);
      best_coordinates[3] = 350*cos(theta);
      best_coordinates[4] = 350*sin(theta);
    } else {
      best_coordinates[3] = best_coordinates[3]*350;
      best_coordinates[4] = best_coordinates[4]*350;
    }
    // More debuggins statements
    //std::cout << "GrDe X: " << best_coordinates[3] << "   GrDe Y: " << best_coordinates[4] << std::endl;
    //std::cout << "TrY: " << eventdata->fGenY/357 << "  muY: " << mu[1] << "  GmuY: " << gradient_mu[1] <<  std::endl;
    //std::cout << "Z: " << eventdata->fGenZ << "  Asym: " << asym << "  Top Hits: " << tb_all[0] << "  Bot Hits: " << tb_all[1] << std::endl;
    //std::cout << "Error: " << error << std::endl;
    //std::cout << "lr_a: " << lr_a << "  lr_sig: " << lr_sigma << "  lr mu x: " << lr_mu_x << "  lr mu y: " << lr_mu_y << std::endl;
    
    
    // Make sure that the reconstructed r is not outside of the detector!!!
    // ---------------------------------------------------------------------------
    // Set best coordinates as reconstructed vertex 
    // ---------------------------------------------------------------------------
    reco_vertex[0] = best_coordinates[3];
    reco_vertex[1] = best_coordinates[4];	

    //--------------------------
    //POST FIRST RUN CORRECTIONS
    //-------------------------
    
    std::cout << "True energy: " << energy << std::endl;

    //------------------------------------------------------------------------------
    //Reconstruct the energy using tsh, the total number of scintillaton hits
    //------------------------------------------------------------------------------
     double tsh = tb_all[0]+tb_all[1];

     double reco_e = 4.84643*(1/pow(10,13))*pow(tsh,3) - 4.8728*(1/pow(10,10))*pow(tsh,2) + 0.00202442*tsh - 0.0569654;

     //figure out which of the 10 simulated energies the reco is closest to and correct accordingly 
     std::vector<double> e_cutoff = {7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5};
     
     bool sorted = false;
     for (int c=0;c<9;c++) {
       if (reco_e<e_cutoff[c] and sorted==false) {
	 energy=5*c+5;
	 sorted=true;
       }
     }
     if (sorted==false) {
       energy=50;
       sorted=true;
     }
     std::string reco_energy = std::to_string(energy);


    //-------------------------------------------------
    //Open text file corresponding to energy and grab correction data
    //------------------------------------------------
    
    std::vector<std::vector<double>> mean_z_error(16, std::vector<double> (16,0));
    std::vector<std::vector<double>> mean_r_error(16, std::vector<double> (16,0));

    std::cout << "Reco energy: " << energy << std::endl;
    std::string line;
    std::string fname = "/cluster/tufts/wongjiradlabnu/bharri06/data/corr_data/data_";
    fname += reco_energy;
    fname += "mev.txt";
    std::ifstream mefile (fname); 
    if ( mefile.is_open()) {
      //skip header
      getline(mefile,line);
      //grab mean z data
      for (int z=0;z<16;z++) {
	getline(mefile,line);
	for (int r=0;r<16;r++) {
	  mean_z_error[r][z]=std::stod( line.substr(0, line.find(' ')));
	  line.erase(0, line.find(' ')+1);
	  //debug
	  //std::cout << mean_z_error[r][z] << " ";
	}
	//std::cout << std::endl;
      }
      std::cout << std::endl;
      //skip blank and second header
      getline(mefile,line);
      getline(mefile,line);
      //grab mean r data
      for (int z=0;z<16;z++) {
        getline(mefile,line);
        for (int r=0;r<16;r++) {
          mean_r_error[r][z]=std::stod( line.substr(0, line.find(' ')));
          line.erase(0, line.find(' ')+1);
          //for debugging
          //std::cout << mean_r_error[r][z] << " ";
	}
	//std::cout << std::endl;
      }
      mefile.close();
    } else {
      std::cout << "unable to open correction file" << std::endl;
    }
    
   //------------------------------
   //Correct z coordinate of vertex
   //-----------------------------
   // loop through all the zones
   for (int r=0; r<350; r+=21.875) {
     int rindex = r/21.875;
     for (int z=-400; z<400; z+=50) {
       int zindex = (z+400)/50;
       // see which are in current zone
       if (reco_vertex[3] > r && reco_vertex[3] < r+21.875) {
         if (reco_vertex[2] > z && reco_vertex[2] < z+50) {
           reco_vertex[2]-=mean_z_error[rindex][zindex];
         }
       }
     } //end of z zone loop
   } //end or r zone loop  

   //------------------------------
   //Correct r coordinate of vertex
   //-----------------------------
    // loop through all the zones
   for (int r=0; r<350; r+=21.875) {
     int rindex = r/21.875;
     for (int z=-400; z<400; z+=50) {
       int zindex = (z+400)/50;
       // see which are in current zone
       if (reco_vertex[3] > r && reco_vertex[3] < r+21.875) {
         if (reco_vertex[2] > z && reco_vertex[2] < z+50) {
           reco_vertex[3]-=mean_r_error[rindex][zindex];
         }
       }
     } //end of z zone loop
   } //end or r zone loop 
    
   // -------------------------------------------------------------------
   // Calculate new x and y coordinates using the corrected r
   // -------------------------------------------------------------------
   //find theta of r 
   double thet=atan2(reco_vertex[1],reco_vertex[0]);
   reco_vertex[0]=reco_vertex[3]*cos(thet);
   reco_vertex[1]=reco_vertex[3]*sin(thet);
 

   //---------------------------------
   //cout and write more stuff to output file 
   //---------------------------------
   
    std::cout << "Corrected Gradient Descent Vertex Error: " << sqrt( pow(reco_vertex[0]-eventdata->fGenX,2) + pow(reco_vertex[1]-eventdata->fGenY,2) + pow(reco_vertex[2]-eventdata->fGenZ,2) ) << std::endl;
    // Write gradient descent reconstructed vertex to outputfile
    outputfile << reco_vertex[0] << " " << reco_vertex[1] << " " << reco_vertex[2] << " ";
  
   //----------------------------
   //DIRECTION RECONSTRUCTION
   //---------------------------
   //Find the direction of the Cherenkov ring        //true = 0, reco = 1
    double highscore[2] = {0, 0};
    double best_x_rot[2] = {0, 0};
    double best_z_rot[2] = {0, 0};
    double best_result[2][3] = {{0, 0, 0}, {0, 0, 0}};
    double xdiff[2] = {0, 0};
    double ydiff[2] = {0, 0};
    double zdiff[2] = {0, 0};
    double HitPosR;
    //increment rotation around x and z axis by 1 degree
    for (double x_rot = 0; x_rot < M_PI; x_rot+=(M_PI/180)) {
        for (double z_rot = 0; z_rot < 2*M_PI; z_rot+=(M_PI/180)) {
            //rotate zenith
            double* temp = rotate_around_x(0, 0, 1, x_rot);
            double* result = rotate_around_z(temp[0], temp[1], temp[2], z_rot);
            //loop through every cherenkov hit, score it, and add it to the total score for this rotation
            double score[2] = {0, 0};
            int size = ceren_hit_x.size();
            for (int i = 0; i < size; i++) {
                //only use sipm ceren hits
                if (channel_tracker[i] > 114) {
                    //find unit vector of hit from true and reconstructed vertex
                    xdiff[0] = {ceren_hit_x[i] - eventdata->fGenX};
     		    xdiff[1] = {ceren_hit_x[i] - reco_vertex[0]};
                    ydiff[0] = {ceren_hit_y[i] - eventdata->fGenY};
                    ydiff[1] = {ceren_hit_y[i] - reco_vertex[1]};
                    zdiff[0] = {ceren_hit_z[i] - eventdata->fGenZ};
                    zdiff[1] = {ceren_hit_z[i] - reco_vertex[2]};
                    HitPosR = sqrt(ceren_hit_x[i]*ceren_hit_x[i]+ceren_hit_y[i]*ceren_hit_y[i]);
                    double p[2] = {sqrt(pow(xdiff[0],2)+pow(ydiff[0],2)+pow(zdiff[0],2)), sqrt(pow(xdiff[1],2)+pow(ydiff[1],2)+pow(zdiff[1],2))};
                    //find theta and phi of hit from new zenith for true and reconstructed vertex
                    double theta[2] = {0, 0};
                    theta[0] = {angle_between(result[0], result[1], result[2], xdiff[0]/p[0], ydiff[0]/p[0], zdiff[0]/p[0])};
                    theta[1] = {angle_between(result[0], result[1], result[2], xdiff[1]/p[1], ydiff[1]/p[1], zdiff[1]/p[1])};
                    double phi[2] = {acos(result[0]/(sin(theta[0]))), acos(result[0]/(sin(theta[1])))};
                    //correct phi when ydiff is negative
                    if (result[1] < 0) {
                        double diff[2] = {M_PI - phi[0], M_PI - phi[1]};
                        phi[0] += 2*diff[0];
                        phi[1] += 2*diff[1];
                    }
		    //give hit a score based on expected theta distribution from TH1D
                    double temp_score[2] = {exp(-pow(theta[0]-0.635,2)/0.1) + exp(-pow(theta[0]-2.3,2)/0.1)/3, exp(-pow(theta[1]-0.635,2)/0.1) + exp(-pow(theta[1]-2.3,2)/0.1)/3};
                    //weight based on accurracy of pmt and sipm
                    if (channel_tracker[i] < 234) {
                        //hit is on an endcap sipm
                        score[0] += 3*temp_score[0];
                        score[1] += 3*temp_score[1];
                        std::cerr << "top sipm hit" << '\n';
                    } else {
                        //hit is on the wall sipm
                        score[0] += temp_score[0];
                        score[1] += temp_score[1];
                    }

                }
                    //using pmt centers, (but not using pmt anymore)
                    // xdiff[0] = {pmt_center[channel_tracker[i]][0] - eventdata->fGenX};
                    // xdiff[1] = {pmt_center[channel_tracker[i]][0] - reco_vertex[0]};
                    // ydiff[0] = {pmt_center[channel_tracker[i]][1] - eventdata->fGenY};
                    // ydiff[1] = {pmt_center[channel_tracker[i]][1] - reco_vertex[1]};
                    // zdiff[0] = {pmt_center[channel_tracker[i]][2] - eventdata->fGenZ};
                    // zdiff[1] = {pmt_center[channel_tracker[i]][2] - reco_vertex[2]};
                    // HitPosR = sqrt(pmt_center[channel_tracker[i]][0]*pmt_center[channel_tracker[i]][0]+pmt_center[channel_tracker[i]][1]*pmt_center[channel_tracker[i]][1]);
            }
            // //fill graphs
            hxzrot_true.Fill( x_rot*180/M_PI, z_rot*180/M_PI, score[0] );
            hxzrot_reco.Fill( x_rot*180/M_PI, z_rot*180/M_PI, score[1] );
            for (int j = 0; j < 2; j++) {
                if (score[j] > highscore[j]) {
                    highscore[j] = score[j];
                    best_x_rot[j] = x_rot;
                    best_z_rot[j] = z_rot;
                    best_result[j][0] = result[0];
                    best_result[j][1] = result[1];
                    best_result[j][2] = result[2];
                }
            }

        }
    }
    // Get the bet reco results
    double xdir = best_result[1][0];
    double ydir = best_result[1][1];
    double zdir = best_result[1][2];

    // Get true directions
    double xtrue = eventdata->fGenXdir;
    double ytrue = eventdata->fGenYdir;
    double ztrue = eventdata->fGenZdir;

    double angle_off = 0;
    //angle_off[0] = {angle_between(best_result[0][0], best_result[0][1], best_result[0][2], xdir, ydir, zdir)};
    angle_off = angle_between(xdir, ydir, zdir, xtrue, ytrue, ztrue);

    std::cout << xdir << " " << ydir << " " << zdir << " " << xtrue << " " << ytrue << " " << ztrue << std::endl;
    std::cout << "Angle off: " << angle_off*180/M_PI << std::endl;

    //print results to outputfile
    outputfile << xtrue << " " << ytrue << " " << ztrue << " " << xdir << " " << ydir << " " << zdir << " " << angle_off*180/M_PI << std::endl;
   


    
    hxzrot_true.Write();
    hxzrot_reco.Write(); 

  } // end of event loop

  outputfile.close();
  
  finput->Close();
  fout->Close();

  return 0;
}


double* rotate_around_x(double x, double y, double z, double rotation) {
    //initialize matricies to store results
    static double result[3] = {0, 0, 0};
    double atemp[3] = {0, 0, 0};
    double btemp[3] = {0, 0, 0};
    //initialize rotation matrix
    double rotation_matrix[3][3] = {{0, 0, 0},
                                    {0, 0, -1},
                                    {0, 1, 0}};
    //create initial vector
    double vector[3] = {x, y, z};
    //rotate around x axis
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            atemp[i] += rotation_matrix[i][k] * vector[k];
        }
    }
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            btemp[i] += rotation_matrix[i][k] * atemp[k];
        }
    }
    for (int i = 0; i < 3; i++) {
        atemp[i] = atemp[i]*sin(rotation);
        btemp[i] = btemp[i]*(1-cos(rotation));
        result[i] = vector[i] + atemp[i] + btemp[i];
    }
    return result;
}


double* rotate_around_y(double x, double y, double z, double rotation) {
    //initialize matricies to store results
    static double result[3] = {0, 0, 0};
    double atemp[3] = {0, 0, 0};
    double btemp[3] = {0, 0, 0};
    //initialize rotation matrix
    double rotation_matrix[3][3] = {{0, 0, 1},
                                    {0, 0, 0},
                                    {-1, 0, 0}};
    //create initial vector
    double vector[3] = {x, y, z};
    //rotate around y axis
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            atemp[i] += rotation_matrix[i][k] * vector[k];
        }
    }
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            btemp[i] += rotation_matrix[i][k] * atemp[k];
        }
    }
    for (int i = 0; i < 3; i++) {
        atemp[i] = atemp[i]*sin(rotation);
        btemp[i] = btemp[i]*(1-cos(rotation));
        result[i] = vector[i] + atemp[i] + btemp[i];
    }
    return result;
}


double* rotate_around_z(double x, double y, double z, double rotation) {
    //initialize matricies to store results
    static double result[3] = {0, 0, 0};
    double atemp[3] = {0, 0, 0};
    double btemp[3] = {0, 0, 0};
    //initialize rotation matrix
    double rotation_matrix[3][3] = {{0, -1, 0},
                                    {1, 0, 0},
                                    {0, 0, 0}};
    //create initial vector
    double vector[3] = {x, y, z};
    //rotate around z axis
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            atemp[i] += rotation_matrix[i][k] * vector[k];
        }
    }
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            btemp[i] += rotation_matrix[i][k] * atemp[k];
        }
    }
    for (int i = 0; i < 3; i++) {
        atemp[i] = atemp[i]*sin(rotation);
        btemp[i] = btemp[i]*(1-cos(rotation));
        result[i] = vector[i] + atemp[i] + btemp[i];
    }
    return result;
}

double angle_between(double x1, double y1, double z1, double x2, double y2, double z2) {
    double dot = x1*x2 + y1*y2 + z1*z2;
    double mag1 = sqrt( pow(x1,2) + pow(y1,2) + pow(z1,2) );
    double mag2 = sqrt( pow(x2,2) + pow(y2,2) + pow(z2,2) );
    return acos(dot/(mag1*mag2));
}
