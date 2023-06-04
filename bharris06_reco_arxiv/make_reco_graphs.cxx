// c++ headers
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>

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

int main( int nargs, char** argv ) {

  std::cout << "Plot Accurracy from data in text file" << std::endl;
  if ( nargs!=2 ) {
    std::cout << "usage: plot_accurracy [outputfile]" << std::endl;
    return 0;
  }

  std::string output_root     = argv[1];

  std::ifstream f(output_root.c_str());

  if ( f.good() ) {
    std::cout << "Output file already exists. Remove first." << std::endl;
    return 0;
  }
  f.close();

  //TApplication app( "app", &nargs, argv );
 
  // ------------------------------------------------------------
  // OUTPUT FILE

  // create root file
  TFile* fout = new TFile( output_root.c_str(), "NEW" );
  
  // declare histograms
  //TCanvas *c1 = new TCanvas("c1","multigraph",500, 500);
  //c1->SetGrid();
  //TMultiGraph *mg = new TMultiGraph();
  

  //std::stringstream hxerror_name;
  //hxerror_name << "hxerror";
  //TH2D hxerror( hxerror_name.str().c_str(), "4 Most Hit PMT; RecoX; RecoX - TrueX", 1000, -400, 400 , 1000, -200, 200);
  
  //std::stringstream hyerror_name;
  //hyerror_name << "hyerror";
  //TH2D hyerror( hyerror_name.str().c_str(), "4 Most Hit PMT; RecoY; RecoY - TrueY", 1000, -400, 400 , 1000, -200, 200);
  
  //std::stringstream hzerror_name;
  //hzerror_name << "hzerror";
  //TH2D hzerror( hzerror_name.str().c_str(), "4 Most Hit PMT; RecoZ; RecoZ - TrueZ", 1000, -400, 400 , 1000, -200, 200);

  // ------------------------------------------------------------
  // INPUT FILE
  
  // variables read in from file
  std::vector<double> energy;
  std::vector<double> reco_vertex[3];
  std::vector<double> true_vertex[3];
  std::vector<double> angle_off;

  // variables created
  std::vector<double> vert_diff[10];
  std::vector<double> r[10];

  std::vector<double> mean_xvert = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<double> sd_xvert = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<double> mean_yvert = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<double> sd_yvert = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<double> mean_zvert = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<double> sd_zvert = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<double> mean_vert = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<double> sd_vert = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<double> mean_angle_off = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<double> sd_angle_off = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
 
  for ( int x=5; x<51; x+=5) {
  	for ( int y=0; y<10; y++) {
		// open text file containing data
		std::string line;
		std::string filename = "/cluster/tufts/wongjiradlab/bharri06/datadir/data_";
		filename += std::to_string(x);
		filename += "_";
		filename += std::to_string(y);
		filename += ".txt";
		//std::cout << filename << std::endl;
		std::ifstream myfile( filename);
		if ( myfile.is_open()) {
			while ( getline( myfile, line)) {
				std::cout << line << std::endl;
				// read in variables from file
				// energy
				energy.push_back( std::stod( line.substr(0, line.find(' '))));
				line.erase(0, line.find(' ')+1);
				// true vertex
				true_vertex[0].push_back( std::stod( line.substr(0, line.find(' '))));
				line.erase(0, line.find(' ')+1);
				true_vertex[1].push_back( std::stod( line.substr(0, line.find(' '))));
                                line.erase(0, line.find(' ')+1);
				true_vertex[2].push_back( std::stod( line.substr(0, line.find(' '))));
                                line.erase(0, line.find(' ')+1);
				// initial reco vertex	




				// z corrected reco vertex
				reco_vertex[0].push_back( std::stod( line.substr(0, line.find(' '))));
                                line.erase(0, line.find(' ')+1);
				reco_vertex[1].push_back( std::stod( line.substr(0, line.find(' '))));
                                line.erase(0, line.find(' ')+1);
				reco_vertex[2].push_back( std::stod( line.substr(0, line.find(' '))));
                        	
				// z and r corrected reco vertex
				
				std::cout << energy.back() << " " << reco_vertex[0].back() << " " << reco_vertex[1].back() << " " << reco_vertex[2].back()
                                          << " "  << true_vertex[0].back() << " " << true_vertex[1].back() << " " << true_vertex[2].back() << " " << angle_off.back() << std::endl;

				// fill error graphs		



				
				// use those to create difference variables which are sorted by energy
				vert_diff[energy.back()].push_back( sqrt( pow(xvert_diff[energy.back()].back(),2) + pow(yvert_diff[energy.back()].back(),2) + pow(zvert_diff[energy.back()].back(),2)));
			 	r.push_back( sqrt( pow(reco_vertex[0].back(),2) + pow(reco_vertex[1].back(),2) + pow(reco_vertex[2].back(),2)));
				//std::cout << energy.back() << " " << xvert_diff.back() << " " << yvert_diff.back() << " " << zvert_diff.back() << " " << vert_diff.back() << std::endl;
			}
			myfile.close();	
		} else { 
			std::cout << "Unable to open file." << std::endl;
		}
	}
  }





 for (int i=0; i<10; i++) {
 	std::cout << "Size - should all be 1000: " << vert_diff[i].size() << std::endl;
 }

 //normal vectors where zone doens't matter
 std::vector<double> vert(10, 0);
 std::vector<int> num_vert(10, 0);
 
 std::vector<double> ms_vert(10, 0);
 std::vector<double> num_ms_vert(10, 0);
 std::vector<double> list_ms_vert[10];
 std::vector<double> sd_ms_vert(10, 0);

 //vectors of vectors where zone matters
 std::vector<std::vector<double>> xvals(10, std::vector<double>(40, 0));
 std::vector<std::vector<double>> xmean(10, std::vector<double>(40, 0));
 std::vector<std::vector<int>> num_xhits_in_zone(10, std::vector<int>(40, 0));

 std::vector<std::vector<double>> yvals(10, std::vector<double>(40, 0));
 std::vector<std::vector<double>> ymean(10, std::vector<double>(40, 0));
 std::vector<std::vector<int>> num_yhits_in_zone(10, std::vector<int>(40, 0));

 std::vector<std::vector<double>> zvals(10, std::vector<double>(40, 0));
 std::vector<std::vector<double>> zmean(10, std::vector<double>(40, 0));
 std::vector<std::vector<int>> num_zhits_in_zone(10, std::vector<int>(40, 0));

 // for (int i=0; i<80; i++) {
//	xvals[i] = 0.0;
//	xmean[i] = 0.0;
 //}
 //
 //
 
 std::cout << "made vectors" << std::endl;
 int energy_index = 0;
 int zone_index = 0;
 //loop through zones
 for (int x=-400; x<401; x+=20) {
	//set xvals
	for (int a=0; a<10; a++) {
		xvals[a][zone_index] = x+10;
		yvals[a][zone_index] = x+10;
		zvals[a][zone_index] = x+10;
	}
	std::cout << "made xvals" << std::endl;
	//std::vector<int> num_hits_in_zone(10, 0);
	//loop through all hits
	for (int i=0; i<10000; i++) {
		

		//assign energy index based on line in file
		for (int num=0; num<10000; num+=1000) {
			if (i>num && i<num+1000) {
				energy_index = (num/1000);
			}
		}

		
		//stuff where zone doesn't matter
		//if (energy_index == 9) {
		//	std::cout << "Reco Vertex: " << reco_vertex[0].at(i) << " " << reco_vertex[1].at(i) << " " << reco_vertex[2].at(i) << std::endl;
		//	std::cout << "True Vertex: " << true_vertex[0].at(i) << " " << true_vertex[1].at(i) << " " << true_vertex[2].at(i) << std::endl;
		//	std::cout << std::endl;
		//}
		//
		// sum vertex differences and keep track of the number of events added to get ready to make the mean
		vert[energy_index] += vert_diff.at(i);
		num_vert[energy_index] += 1;

		// sum medium sphere vertex differences and keep track of the number of events added to get ready to make the mean
		if ( r.at(i)<250) {
			std::cout << "Energy Index: " << energy_index << std::endl;
                	ms_vert[energy_index] += vert_diff.at(i);
                	num_ms_vert[energy_index] += 1;
			std::cout << "Numer sor far: " << num_ms_vert[energy_index] << std::endl; 
			
			//sd_ms_vert[energy_index] += pow(
			list_ms_vert[energy_index].push_back( vert_diff.at(i) );
		}


		

		




		//stuff where zone matters
		//if hit is in zone add it to mean for that energy index and keep track of number of hits
		if ( reco_vertex[0].at(i) > x && reco_vertex[0].at(i) < x+20) {
			xmean[energy_index][zone_index] = xmean[energy_index][zone_index] + (reco_vertex[0].at(i)-true_vertex[0].at(i));
			num_xhits_in_zone[energy_index][zone_index] += 1;
			//std::cout << num_hits_in_zone[energy_index][zone_index] << " " << xmean[energy_index][zone_index] << " " << energy_index << " " << zone_index << std::endl;
		}
		if ( reco_vertex[1].at(i) > x && reco_vertex[1].at(i) < x+20) {
                        ymean[energy_index][zone_index] = ymean[energy_index][zone_index] + (reco_vertex[1].at(i)-true_vertex[1].at(i));
                        num_yhits_in_zone[energy_index][zone_index] += 1;
                        //std::cout << num_hits_in_zone[energy_index][zone_index] << " " << xmean[energy_index][zone_index] << " " << energy_index << " " << zone_index << std::endl;
                }
		if ( reco_vertex[2].at(i) > x && reco_vertex[2].at(i) < x+20) {
                        zmean[energy_index][zone_index] = zmean[energy_index][zone_index] + (reco_vertex[2].at(i)-true_vertex[2].at(i));
                        num_zhits_in_zone[energy_index][zone_index] += 1;
                        //std::cout << num_hits_in_zone[energy_index][zone_index] << " " << xmean[energy_index][zone_index] << " " << energy_index << " " << zone_index << std::endl;
                }
	
	
	}
	//if (num_hits_in_zone[energy_index][zone_index] != 0) {
	//	xmean[energy_index][zone_index] = xmean[energy_index][zone_index]/num_hits_in_zone[energy_index][zone_index];
	//}
	zone_index++;
 }
  // compute all means by dividing by the number of entries 
  for (int i=0; i<10; i++) {
	vert[i] = vert[i]/num_vert[i];
	ms_vert[i] = ms_vert[i]/num_ms_vert[i];

	std::cout << "MS num events kept: " << num_ms_vert[i] << std::endl;
	for (int j=0; j<40; j++) {
		if (num_xhits_in_zone[i][j] != 0) {
                	xmean[i][j] = xmean[i][j]/num_xhits_in_zone[i][j];
		}
		if (num_yhits_in_zone[i][j] != 0) {
                        ymean[i][j] = ymean[i][j]/num_yhits_in_zone[i][j];
                }
		if (num_zhits_in_zone[i][j] != 0) {
                        zmean[i][j] = zmean[i][j]/num_zhits_in_zone[i][j];
                }
	}
  }
  // calculate standard deviation
  for (int energy_index=0; energy_index<10; energy_index++) {
	std::cout << "Standard Deviation Error Check --- Size of array going through: " << list_ms_vert[energy_index].size() << std::endl;
//  	for (int i=0; i<10000; i++) {

		//assign energy index based on line in file
   //             	for (int num=0; num<10000; num+=1000) {
 //                       	if (i>num && i<num+1000) {
     //                           	energy_index = (num/1000);
       //         		}
        //		}
	//	if ( r.at(i)<250) {
//			sd_ms_vert[energy_index] = sd_ms_vert[energy_index] + pow( vert_diff.at(i) - ms_vert[energy_index], 2);
        //	}
//	}
//	sd_ms_vert[energy_index] = sqrt( sd_ms_vert[energy_index]/num_ms_vert[energy_index] );




  }

  std::vector<double> sd_e(10, 0); 


  std::cout << "before erase" << std::endl;
  std::vector<double> xmean_final[10];
  std::vector<double> xvals_final[10];
  std::vector<double> ymean_final[10];
  std::vector<double> yvals_final[10];
  std::vector<double> zmean_final[10];
  std::vector<double> zvals_final[10];
  for (int a=0; a<10; a++) {
  	for (int i=0; i<40; i++) {
		if (xmean[a][i] != 0) {
			xmean_final[a].push_back(xmean[a][i]);
			xvals_final[a].push_back(xvals[a][i]);
		}
		if (ymean[a][i] != 0) {
                        ymean_final[a].push_back(ymean[a][i]);
                        yvals_final[a].push_back(yvals[a][i]);
                }
		if (zmean[a][i] != 0) {
                        zmean_final[a].push_back(zmean[a][i]);
                        zvals_final[a].push_back(zvals[a][i]);
                }
  	}
  }
  std::cout << "after erase" << std::endl;

  //for (int a=0; a<10;a++) { 
//	std::cout << a << std::endl;
  //	for (int i=0; i<xmean_final[a].size(); i++) {
//		std::cout << xvals_final[a][i] << " " << xmean_final[a][i] << std::endl;
  //	}
//	std::cout << std::endl;
//	std::cout << std::endl;
  //}
  
  for(int i=0; i<num_ms_vert.size(); i++) {
	num_ms_vert[i] = num_ms_vert[i]/40000;	  
  }



  std::vector<double> evals = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50};

  TCanvas *c6 = new TCanvas("c6", "events_kept", 50, 1);
  c6->SetGrid();
  TGraphErrors *ek = new TGraphErrors(10, evals.data(), num_ms_vert.data()); 
  ek->SetMarkerColor(kBlue);
  ek->SetMarkerStyle(21);
  ek->SetName("ek");
  ek->Draw("apl");
  ek->Write();



  TCanvas *c4 = new TCanvas("c4", "hverror", 50, 60);
  c4->SetGrid();
  TGraphErrors *hverror = new TGraphErrors(10, evals.data(), vert.data());
  
  hverror->SetMarkerColor(kBlue);
  hverror->SetMarkerStyle(21);
  hverror->SetName("hverror");
  hverror->Draw("apl");
  hverror->Write();

  TCanvas *c5 = new TCanvas("c5", "hvmserror", 50, 60);
  c5->SetGrid();
  TGraphErrors *hmsverror = new TGraphErrors(10, evals.data(), ms_vert.data(), sd_e.data(), sd_ms_vert.data());

  hmsverror->SetMarkerColor(kGreen);
  hmsverror->SetMarkerStyle(21);
  hmsverror->SetName("hmsverror");
  hmsverror->Draw("apl");
  hmsverror->Write();

  for (int a=0; a<10; a++) {
	
  std::cout << "0" << std::endl;
  
  std::string add  = std::to_string(a);
  
  //xerror
  std::string xname = "c1";
  std::string xtitle = "hxerror";
  xname += add;
  xtitle += add;

  //yerror
  std::string yname = "c2";
  std::string ytitle = "hyerror";
  yname += add;
  ytitle += add;

  //zerror
  std::string zname = "c3";
  std::string ztitle = "hzerror";
  zname += add;
  ztitle += add;

  TCanvas *c1 = new TCanvas(xname.c_str(), xtitle.c_str(), 800, 400);
  c1->SetGrid();
  TGraphErrors *hxerror = new TGraphErrors(xvals_final[a].size(), xvals_final[a].data(), xmean_final[a].data());
  
  TCanvas *c2 = new TCanvas(yname.c_str(), ytitle.c_str(), 800, 400);
  c2->SetGrid();
  TGraphErrors *hyerror = new TGraphErrors(yvals_final[a].size(), yvals_final[a].data(), ymean_final[a].data());

  TCanvas *c3 = new TCanvas(zname.c_str(), ztitle.c_str(), 800, 400);
  c3->SetGrid();
  TGraphErrors *hzerror = new TGraphErrors(zvals_final[a].size(), zvals_final[a].data(), zmean_final[a].data());


  // Fit Code
    double xmin = xvals[a][0];
    double xmax = xvals[a][xvals.size()-1];
    TF1 myxFunc_fifth("Polynomial5thDegree","pol5",xmin,xmax);
    hxerror->Fit(&myxFunc_fifth);
    double x5 = myxFunc_fifth.GetParameters()[0];
    double x4 = myxFunc_fifth.GetParameters()[1];
    double x3 = myxFunc_fifth.GetParameters()[2];
    double x2 = myxFunc_fifth.GetParameters()[3];
    double x1 = myxFunc_fifth.GetParameters()[4];
    double xc  = myxFunc_fifth.GetParameters()[5];
    TGraph xline(0);
    int xi =0;
    for (int x_coord = xmin; x_coord <= xmax; x_coord++){
      double y= x5*pow(x_coord,5) + x4*pow(x_coord,4) + x3*pow(x_coord,3) + x2*pow(x_coord,2) + x1*x_coord + xc;
      xline.SetPoint(xi,x_coord,y);
      xi++;
    }

    double ymin = yvals[a][0];
    double ymax = yvals[a][yvals.size()-1];
    TF1 myyFunc_fifth("Polynomial5thDegree","pol5",ymin,ymax);
    hyerror->Fit(&myyFunc_fifth);
    double y5 = myyFunc_fifth.GetParameters()[0];
    double y4 = myyFunc_fifth.GetParameters()[1];
    double y3 = myyFunc_fifth.GetParameters()[2];
    double y2 = myyFunc_fifth.GetParameters()[3];
    double y1 = myyFunc_fifth.GetParameters()[4];
    double yc  = myyFunc_fifth.GetParameters()[5];
    TGraph yline(0);
    int yi =0;
    for (int y_coord = ymin; y_coord <= ymax; y_coord++){
      double y= y5*pow(y_coord,5) + y4*pow(y_coord,4) + y3*pow(y_coord,3) + y2*pow(y_coord,2) + y1*y_coord + yc;
      yline.SetPoint(yi,y_coord,y);
      yi++;
    }

    double zmin = zvals[a][0];
    double zmax = zvals[a][zvals.size()-1];
    TF1 myzFunc_fifth("Polynomial5thDegree","pol5",zmin,zmax);
    hzerror->Fit(&myzFunc_fifth);
    double z5 = myzFunc_fifth.GetParameters()[0];
    double z4 = myzFunc_fifth.GetParameters()[1];
    double z3 = myzFunc_fifth.GetParameters()[2];
    double z2 = myzFunc_fifth.GetParameters()[3];
    double z1 = myzFunc_fifth.GetParameters()[4];
    double zc  = myzFunc_fifth.GetParameters()[5];
    TGraph zline(0);
    int zi =0;
    for (int z_coord = zmin; z_coord <= zmax; z_coord++){
      double y= z5*pow(z_coord,5) + z4*pow(z_coord,4) + z3*pow(z_coord,3) + z2*pow(z_coord,2) + z1*z_coord + zc;
      zline.SetPoint(zi,z_coord,y);
      zi++;
    }
  // End Fit Code

  std::cout << "3" << std::endl;
  hxerror->SetMarkerColor(kBlue);
  hxerror->SetMarkerStyle(21);

  hyerror->SetMarkerColor(kRed);
  hxerror->SetMarkerStyle(21);

  hzerror->SetMarkerColor(kGreen);
  hzerror->SetMarkerStyle(21);

  if (a==0) {
	hxerror->SetName("hxerror5mev");
	hyerror->SetName("hyerror5mev");
	hzerror->SetName("hzerror5mev");
  } else if (a==1) {
	hxerror->SetName("hxerror10mev");
	hyerror->SetName("hyerror10mev");
	hzerror->SetName("hzerror10mev");
  } else if (a==2) {
	hxerror->SetName("hxerror15mev");
	hyerror->SetName("hyerror15mev");
	hzerror->SetName("hzerror15mev");
  } else if (a==3) {
	hxerror->SetName("hxerror20mev");
	hyerror->SetName("hyerror20mev");
	hzerror->SetName("hzerror20mev");
  } else if (a==4) {
        hxerror->SetName("hxerror25mev");
	hyerror->SetName("hyerror25mev");
	hzerror->SetName("hzerror25mev");
  } else if (a==5) {
        hxerror->SetName("hxerror30mev");
	hyerror->SetName("hyerror30mev");
	hzerror->SetName("hzerror30mev");
  } else if (a==6) {
        hxerror->SetName("hxerror35mev");
	hyerror->SetName("hyerror35mev");
	hzerror->SetName("hzerror35mev");
  } else if (a==7) {
        hxerror->SetName("hxerror40mev");
	hyerror->SetName("hyerror40mev");
	hzerror->SetName("hzerror40mev");
  } else if (a==8) {
        hxerror->SetName("hxerror45mev");
	hyerror->SetName("hyerror45mev");
	hzerror->SetName("hzerror45mev");
  } else if (a==9) {
        hxerror->SetName("hxerror50mev");
	hyerror->SetName("hyerror50mev");
	hzerror->SetName("hzerror50mev");
  }
   hxerror->Draw("apl");
   hxerror->Write(); 
   hyerror->Draw("apl");
   hyerror->Write();
   hzerror->Draw("apl");
   hzerror->Write();
  }


  std::cout << "Made it out" << std::endl;
   
  //hxerror.Write();
  //hyerror.Write();
  //hzerror.Write();

  // make histograms
  //for (int i=0; i<mean_xvert.size(); i++) {
  //	std::cout << mean_angle_off.at(i) << " ";
  //}



  fout->Close();
  
  return 0;
}
