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
#include "TFile.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TAttLine.h"
#include "Rtypes.h"

// ----------------------------------------------------------------
// -------- VISUALIZING CHAOS WITH THE BUTTERFLY EFFECT -----------
// ----------------------------------------------------------------
//  Ben Harris - May 10, 2022 - Computational Physics Final Project
// ----------------------------------------------------------------

int main( int nargs, char** argv ) {
  
  // --------------
  // ARGUMENT CHECK
  // --------------
  
  std::cout << "Visualizing Chaos..." << std::endl;
  if ( nargs!=2 ) {
    std::cout << "usage: chaos [outputfile]" << std::endl;
    return 0;
  }
  // get output rootfile as a string
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
 
  // ------------------------------------------------
  // MODEL LORENZ ATTRACTOR WITH CONVECTION EQUATIONS
  // ------------------------------------------------
  
  //constant parameters that determine the properties of the atmosphere
  double a = 10;
  double b = 8/3;
  double c = 28;

  //set timestep
  double dt = 0.0001; 				//0.0001 is ideal

  //set desired number of iterations
  const int iterations = 500000;  		//1,000,000 at a dt of 0.0001 is ideal
  int index = 0;

  //initialize vectors to store results for 3 different sets of initial conditions
  std::vector<std::vector<double>> x (3, std::vector<double> (iterations,0));
  std::vector<std::vector<double>> y (3, std::vector<double> (iterations,0));
  std::vector<std::vector<double>> z (3, std::vector<double> (iterations,0));
 
  //set initial coordinates for blue 
  x[0][0] = 5;
  y[0][0] = 7;
  z[0][0] = 10;

  //set initial coordinates for green
  x[1][0] = 5.000001;
  y[1][0] = 7.000001;
  z[1][0] = 10.000001;

  //set initial coordinates for red
  x[2][0] = 4.999999;
  y[2][0] = 6.999999;
  z[2][0] = 9.999999;

  //make canvas
  TCanvas *c1 = new TCanvas("c1","Lorenz Attractor",40,50);
  
  //debugging
  std::cout << "Iterations: " << iterations << "   Time step: " << dt << "   Max simulation time: " << iterations*dt << std::endl;

  //loop that updates x,y,z for the desired number of iterations
  while (index<iterations-1) {
    //update position for each color
    for (int color=0;color<3;color++) {
      x[color][index+1] = x[color][index] - a*x[color][index]*dt + a*y[color][index]*dt;
      y[color][index+1] = y[color][index] - x[color][index]*z[color][index]*dt + c*x[color][index]*dt - y[color][index]*dt;
      z[color][index+1] = z[color][index] + x[color][index]*y[color][index]*dt - b*z[color][index]*dt;
    }
    //go to next iteration
    index++;
  }
 
  //more debugging
  std::cout << "Simulation stopped after " << index+1 << "/" << iterations << " iterations" << std::endl;
  
  //cut off the last 30k iterations for each color
  //used (not much) for debugging and visualizing certain sections
  //for (int color=0;color<3;color++) {
  //  x[color].erase(x[color].begin(),x[color].end()-30000);
  //  y[color].erase(y[color].begin(),y[color].end()-30000);
  //  z[color].erase(z[color].begin(),z[color].end()-30000);
  //}


  //plot red results
  TGraph* red = new TGraph(iterations,x[2].data(),z[2].data());
  red->SetTitle("Lorenz Attractor; X; Z");
  red->SetName("red");
  red->SetMarkerColor(2);
  red->SetLineWidth(1);
  red->Draw();
  red->Write();

  //plot green results
  TGraph* green = new TGraph(iterations,x[1].data(),z[1].data());
  green->SetTitle("Lorenz Attractor; X; Z");
  green->SetName("green");
  green->SetMarkerColor(3);
  green->SetLineWidth(1);
  green->Draw();
  green->Write();

  //plot blue results
  TGraph* blue = new TGraph(iterations,x[0].data(),z[0].data());
  blue->SetTitle("Lorenz Attractor; X; Z");
  blue->SetName("blue");
  blue->SetMarkerColor(4);
  blue->SetLineWidth(1);
  blue->Draw();
  blue->Write();

  //plot all results
  TMultiGraph* rgb = new TMultiGraph();
  rgb->Add(red);
  rgb->Add(green);
  rgb->Add(blue);
  rgb->SetTitle("The Butterfly Effect; X; Z");
  rgb->SetName("rgb");
  rgb->Draw();
  rgb->Write();

  //zoom in on blue results
  TGraph* zoom = new TGraph(iterations,x[0].data(),z[0].data());
  zoom->SetTitle("Lorenz Attractor; X; Z");
  zoom->SetName("zoom");
  zoom->SetMarkerColor(4);
  zoom->SetLineWidth(1);
  zoom->Draw();
  zoom->GetXaxis()->SetRangeUser(0,15);
  zoom->GetYaxis()->SetRangeUser(10,40);
  zoom->Write();

  //zoom in further on blue results
  TGraph* zoomzoom = new TGraph(iterations,x[0].data(),z[0].data());
  zoomzoom->SetTitle("Lorenz Attractor; X; Z");
  zoomzoom->SetName("zoomzoom");
  zoomzoom->SetMarkerColor(4);
  zoomzoom->SetLineWidth(1);
  zoomzoom->Draw();
  zoomzoom->GetXaxis()->SetRangeUser(6,10);
  zoomzoom->GetYaxis()->SetRangeUser(15,25);
  zoomzoom->Write();
  
  fout->Close();

  return 0;
}
