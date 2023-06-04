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
//#include "TApplication.h"
#include "TMath.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TNamed.h"
#include "TObject.h"
#include "TGraph.h"
#include "TGraph2DErrors.h"
#include "TF1.h"
#include "TGraphErrors.h"

int main( int nargs, char** argv ) {

  std::cout << "Plot Accurracy from data in text file" << std::endl;
  if ( nargs!=2 ) {
    std::cout << "usage: plot_accurracy [outputfile]" << std::endl;
    return 0;
  }

//  TApplication app( "app", &nargs, argv );
 
  // ------------------------------------------------------------
  // OUTPUT FILE

  std::string output_root     = argv[1];
  
  TFile* fout = new TFile( output_root.c_str(), "NEW" );
  
  // declare histograms
  //TCanvas *c1 = new TCanvas("c1","multigraph",500, 500);
  //c1->SetGrid();
  //TMultiGraph *mg = new TMultiGraph();
  
  // ------------------------------------------------------------
  // INPUT FILE
  
  // variables read in from file
  std::vector<double> true_vertex[10][4];
  std::vector<double> reco_vertex[10][4];

  for ( int x=5; x<51; x+=5) {
  	for ( int y=0; y<10; y++) {
		// open text file containing data
		std::string line;
		std::string filename = "/cluster/tufts/wongjiradlabnu/bharri06/data/reco_data/T4data10x/data_";
		filename += std::to_string(x);
		filename += "_";
		filename += std::to_string(y);
		filename += ".txt";
		std::cout << filename << std::endl;
		std::ifstream myfile( filename);
		if ( myfile.is_open()) {
			while ( getline( myfile, line)) {
				// read in variables from file
				int energy = std::stoi( line.substr(0, line.find(' ')));
				int i = (energy-5)/5;
				line.erase(0, line.find(' ')+1);
				true_vertex[i][0].push_back( std::stod( line.substr(0, line.find(' '))));
				line.erase(0, line.find(' ')+1);
				true_vertex[i][1].push_back( std::stod( line.substr(0, line.find(' '))));
				line.erase(0, line.find(' ')+1);
				true_vertex[i][2].push_back( std::stod( line.substr(0, line.find(' '))));
				true_vertex[i][3].push_back(sqrt(pow(true_vertex[i][0].back(),2)+pow(true_vertex[i][1].back(),2)));
				line.erase(0, line.find(' ')+1);
				reco_vertex[i][0].push_back( std::stod( line.substr(0, line.find(' '))));
                                line.erase(0, line.find(' ')+1);
                                reco_vertex[i][1].push_back( std::stod( line.substr(0, line.find(' '))));
                                line.erase(0, line.find(' ')+1);
                                reco_vertex[i][2].push_back( std::stod( line.substr(0, line.find(' '))));
				reco_vertex[i][3].push_back(sqrt(pow(reco_vertex[i][0].back(),2)+pow(reco_vertex[i][1].back(),2)));
			}
			myfile.close();	
		} else { 
			std::cout << "Unable to open file." << std::endl;
		}
	}
  }

  //-----------------------------------------------------------------------------------------------------------------------------------------------------
  //Part negative one
  //Making histogram of reco error as a fucntion of energy
  
  //find mean reco error at each energy
  std::vector<double> initial_error[10];
  std::vector<double> mean_initial_error(10,0);
  for (int i=0; i<10; i++) {
  	int size = true_vertex[i][0].size();
	for (int j=0; j<size; j++) {
		double xdiff = abs(reco_vertex[i][0].at(j)-true_vertex[i][0].at(j));
		double ydiff = abs(reco_vertex[i][1].at(j)-true_vertex[i][1].at(j));
		double zdiff = abs(reco_vertex[i][2].at(j)-true_vertex[i][2].at(j));
		initial_error[i].push_back(sqrt(pow(xdiff,2)+pow(ydiff,2)+pow(zdiff,2)));
		mean_initial_error[i] += initial_error[i].back();
	}
	mean_initial_error[i] = mean_initial_error[i]/size;
  }

  //find standard deviation at each energy
  std::vector<double> blank(10,0);
  std::vector<double> initial_sd(10,0);
  for (int i=0; i<10; i++) {
        int size = initial_error[i].size();
        for (int j=0; j<size; j++) {
                initial_sd[i] += pow(abs(initial_error[i].at(j)-mean_initial_error[i]),2);
        }
        initial_sd[i] = sqrt(initial_sd[i]/size);
  }

  //make energies
  std::vector<double> energies = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50};

  //make graph
  TGraphErrors *c1 = new TGraphErrors( 10, energies.data(), mean_initial_error.data(), blank.data(), initial_sd.data() );
  c1->SetMarkerColor(kBlue);
  c1->SetMarkerStyle(21);
  c1->SetName("initial_reco_rertex_error");
  c1->GetXaxis()->SetTitle("Energies (MeV)");
  c1->GetYaxis()->SetTitle("Mean |Reco V - True V| Error");
  c1->Draw("apl");
  c1->Write();

  //-----------------------------------------------------------------------------------------------------------------------------------------------------
  //Part two
  //Making histogram of reco z error as a function of reco z and reco r

  std::string axis = ";Reco R;Reco Z";
  for (int i=0; i<10; i++) {
        //make histogram header and axis titles
        std::string energy = std::to_string(5*i+5);
        std::string title = energy + "mev" + axis;
        //make histrogram name
        std::stringstream hzerror_name;
        hzerror_name << "hzerror" << energy << "mev";
        //make histogram
        TH2D hzerror( hzerror_name.str().c_str(), title.c_str(), 96, 0, 480, 160, -400, 400);
        //fill histogram
        int size = reco_vertex[i][2].size();
        for (int j=0; j<size; j++) {
                hzerror.Fill(reco_vertex[i][3].at(j), reco_vertex[i][2].at(j), reco_vertex[i][3].at(j)-true_vertex[i][3].at(j));
        }
	hzerror.GetZaxis()->SetRangeUser(-200,200);
	hzerror.Draw("colz");
        hzerror.Write();
  }


  //-----------------------------------------------------------------------------------------------------------------------------------------------------
  //Part three
  //Making plot of mean reco z error as a function of reco z and reco r
 
  std::vector<std::vector<std::vector<double>>> mean_zerror(10, std::vector<std::vector<double>> (20, std::vector<double> (10, 0)));
  std::vector<std::vector<std::vector<double>>> num_in_zone(10, std::vector<std::vector<double>> (20, std::vector<double> (10, 0)));
  std::vector<double> list_z_error[10][20][10];
  std::vector<double> list_reco_z[10][20][10];
  std::vector<double> list_true_z[10][20][10];
  int zbins = 20;
  int zwidth = 800/zbins;
  int rbins = 10;
  int rwidth = 400/rbins;
  //find mean error for each bin
  for (int i=0; i<10; i++) {
  	for (int z=0; z<zbins; z++) {
		for (int r=0; r<rbins; r++) {
			int size = true_vertex[i][0].size();
			for (int j=0; j<size; j++) {
				if (reco_vertex[i][2].at(j) > -400+(zwidth*z) and reco_vertex[i][2].at(j) < -400+(zwidth*(z+1)) ) {
					if(reco_vertex[i][3].at(j) > (rwidth*r) and reco_vertex[i][3].at(j) < rwidth*(r+1) ) {	
						mean_zerror[i][z][r] += reco_vertex[i][2].at(j)-true_vertex[i][2].at(j);
						num_in_zone[i][z][r]++;
						list_z_error[i][z][r].push_back(reco_vertex[i][2].at(j)-true_vertex[i][2].at(j));
						list_reco_z[i][z][r].push_back(reco_vertex[i][2].at(j));
						list_true_z[i][z][r].push_back(true_vertex[i][2].at(j));
					}
				}
			}
			mean_zerror[i][z][r] = mean_zerror[i][z][r]/num_in_zone[i][z][r];
		//	std::cout << "E: " << i << "  Z: " << z << "  R: " << r << "  M: " << mean_zerror[i][z][r] << std::endl;
		}
	}	
  }
  //find sd for each bin
  std::vector<std::vector<std::vector<double>>> sd_z(10, std::vector<std::vector<double>> (20, std::vector<double> (10, 0)));
  for (int i=0; i<10; i++) {
  	for (int z=0; z<zbins; z++) {
		for (int r; r<rbins; r++) {
			int size = list_reco_z[i][z][r].size();
			for (int j=0; j<size; j++) {
				sd_z[i][z][r] += pow(list_reco_z[i][z][r].at(j)-mean_zerror[i][z][r],2);	
			}
			sd_z[i][z][r] = sqrt(sd_z[i][z][r]/size);
		}
	}
  }

  //make recozvals for x and y axis
  std::vector<std::vector<double>> zvals(10, std::vector<double> (20, 0));
  std::vector<std::vector<double>> rvals(10, std::vector<double> (10, 0));
  for (int i=0; i<10; i++) {
  	int rindex = 0;
	for (int r=0; r<400; r+=rwidth) {
        	rvals[i][rindex] = r+(rwidth/2);
       		rindex++;
  	}
	int zindex = 0;
	for (int z=-400; z<400; z+=zwidth) {
		zvals[i][zindex] = z+(zwidth/2);
		zindex++;
	
	}
  }

  //make all nan in  mean z error into zeros to make sure its good for a correction
  for (int i=0; i<10; i++) {
        for (int z=0; z<zbins; z++) {
		for (int r=0; r<rbins; r++) {
        	        if (std::isnan(mean_zerror[i][z][r])) {
               			mean_zerror[i][z][r] = 0;
               		}
        		if (std::isnan(sd_z[i][z][r])) {
				sd_z[i][z][r] = 0;
			}
			std::cout << "E: " << i << "  Z: " << z << "  R: " << r << "  M: " << mean_zerror[i][z][r] << std::endl;
		}
	}
  }


// //erase the blanks by creating final_zvals, final_mean, final_sd
// std::vector<double> xa_zr[10];
// std::vector<double> ya_zr[10];
// std::vector<double> xsd_zr[10];
// std::vector<double> ysd_zr[10];
// //better reco 
// std::vector<double> xa_za[10];
// std::vector<double> ya_za[10];
// std::vector<double> ysd_za[10];
//
// for (int i=0; i<10; i++) {
//        for (int j=0; j<bins; j++) {
//                if ( meanzrerror[i][j] > -500 and meanzrerror[i][j] < 500 and meanzrerror[i][j] != 0) {
//                        xa_zr[i].push_back( recozvals[i][j] );
//                        ya_zr[i].push_back( meanzrerror[i][j] );
//                        xsd_zr[i].push_back( 0 );
//                        ysd_zr[i].push_back( sd_zr[i][j] );
//                }
//        	if( meanzaerror[i][j] > -500 and meanzaerror[i][j] < 500 and meanzaerror[i][j] !=0) {
//			xa_za[i].push_back( recozvals[i][j] );
//			ya_za[i].push_back( meanzaerror[i][j] );
//			xsd_za[i].push_back( 0 );
//			ysd_za[i].push_back( sd_za[i][j] );
//		}
//	}
// } 
 std::vector<std::vector<double>> rblank (10, std::vector<double> (10, 0)); 
 std::vector<std::vector<double>> zblank (10, std::vector<double> (20, 0));
 //make plots
 
 std::string ax = ";Reco R;Reco Z";
  for (int i=0; i<10; i++) {
        //make histogram header and axis titles
        std::string energy = std::to_string(5*i+5);
        std::string title = energy + "mev" + ax;
        //make histrogram name
        std::stringstream hzmeancorr_name;
        hzmeancorr_name << "hzmeancorr" << energy << "mev";
        //make histogram
        TH2D hzmeancorr( hzmeancorr_name.str().c_str(), title.c_str(), 12, 0, 480, 20, -400, 400);
        //fill histogram
	//
	for (int z=0; z<zbins; z++) {
		for (int r=0; r<rbins; r++) {
			hzmeancorr.Fill(rvals[i][r], zvals[i][z], mean_zerror[i][z][r]);
		}
	}
        hzmeancorr.GetZaxis()->SetRangeUser(-100,100);
        hzmeancorr.Draw("colz");
        hzmeancorr.Write();
  }
 
 
 
// for (int a=0; a<10; a++) {
//
//        std::string title  = std::to_string(a*5+5);
//	title += "mev";
//        std::string name = "z_correction";
//	name += title;
//
//
// 	TCanvas *c69 = new TCanvas(name.c_str(), title.c_str(), 500, 800);
// 	c69->SetGrid();
// 	TGraph2DErrors *zbaseline = new TGraph2DErrors(zbins*rbins);
//	for (int z=0; z<zbins; z++) {
//		for (int r=0; r<rbins; r++) {
//			zbaseline->SetPoint(20*z+r, rvals[a][r], zvals[a][z], mean_zerror[a][z][r]);
//			zbaseline->SetPointError(20*z+r, rblank[a][r], zblank[a][z], sd_z[a][z][r]); 
//		}
//	}
//	
//	
//	//, zvals[a].data(), rvals[a].data(), mean_zerror[a].data()); //, blank[a].data(), sd_zr[a].data() );
//	//TGraphErrors *zbaseline = new TGraphErrors(1000, recoz[a].data(), zerror[a].data() );
//
//
//       // // Fit Code
//       // double min = xaxis[a].front();
//       // double max = xaxis[a].back();
//       // TF1 myxFunc_fourth("Polynomial4thDegree","pol4",min,max);
//       // zbaseline->Fit(&myxFunc_fourth);
//       // double x4 = myxFunc_fourth.GetParameters()[1];
//       // double x3 = myxFunc_fourth.GetParameters()[2];
//       // double x2 = myxFunc_fourth.GetParameters()[3];
//       // double x1 = myxFunc_fourth.GetParameters()[4];
//       // double xc  = myxFunc_fourth.GetParameters()[5];
//       // TGraph line(0);
//       // int b=0;
//       // for (int coord = min; coord <= max; coord++){
//       //         double y = x4*pow(coord,4) + x3*pow(coord,3) + x2*pow(coord,2) + x1*coord + xc;
//       //         line.SetPoint(b,coord,y);
//       //         b++;
//       // }
//        // End Fit Code
//
//	zbaseline->SetName(name.c_str());
//	zbaseline->SetTitle(title.c_str());	
//	zbaseline->SetMarkerColor(kBlue);
//        zbaseline->SetMarkerStyle(21);
//        zbaseline->GetXaxis()->SetTitle("Reco Z");
//        zbaseline->GetYaxis()->SetTitle("Reco Z - True Z");
//	//zbaseline->GetYaxis()->SetRangeUser(-100,100);
//	zbaseline->Draw("colz");
//	zbaseline->Write();

	// // Fit Code
        //double bmin = xa[a].front();
        //double bmax = xa[a].back();
        //TF1 mybFunc_fourth("Polynomial4thDegree","pol4",bmin,bmax);
        //bzbaseline->Fit(&mybFunc_fourth);
        //double b4 = mybFunc_fourth.GetParameters()[1];
        //double b3 = mybFunc_fourth.GetParameters()[2];
        //double b2 = mybFunc_fourth.GetParameters()[3];
        //double b1 = mybFunc_fourth.GetParameters()[4];
        //double bc  = mybFunc_fourth.GetParameters()[5];
        //TGraph bline(0);
        //int p=0;
        //for (int coord = bmin; coord <= bmax; coord++){
        //        double c = b4*pow(coord,4) + b3*pow(coord,3) + b2*pow(coord,2) + b1*coord + bc;
        //        bline.SetPoint(p,coord,c);
        //        p++;
        //}
        //// End Fit Code
	//
 //}
 


  //-----------------------------------------------------------------------------------------------------------------------------------------------------
  //Part four
  //Correct reco z and make final graph of reco z error as function of energy
  

  //correct data
  for (int i=0; i<10; i++) {
  	for (int z=0; z<zbins;z++) {
		for (int r=0; r<rbins;r++) {
			int size = reco_vertex[i][2].size();
			for (int j=0; j<size; j++) {
				if (reco_vertex[i][2].at(j) > -400+(zwidth*z) and reco_vertex[i][2].at(j) < -400+(zwidth*(z+1)) ) {
                                        if(reco_vertex[i][3].at(j) > (rwidth*r) and reco_vertex[i][3].at(j) < rwidth*(r+1) ) {
                                		reco_vertex[i][2].at(j) -= mean_zerror[i][z][r];
					}
                                }
			}
		}
	}
  }

  //first correct reco z data
  
  // for (int i =0; i<10; i++) {
 // 	if (i==0) {
 //       	for (int j=0; j<1000; j++) {
 //       		double correction = -21.6288 + 0.0304921*recoz[i].at(j) + 0.000848095*pow(recoz[i].at(j),2) - 1.40038*(1/pow(10,6))*pow(recoz[i].at(j),3) - 2.94491*(1/pow(10,9))*pow(recoz[i].at(j),4);
 //       		recoz[i].at(j) -= correction;
 //       		double all = -21.9561 + 0.0304921*recoz[i].at(j) + 0.000821394*pow(recoz[i].at(j),2) -1.58377*(1/pow(10,6))*pow(recoz[i].at(j),3) -2.42671*(1/pow(10,9))*pow(recoz[i].at(j),4);                    
 //       		bz[i].at(j) = recoz[i].at(j) - all;
 //       	}
 //       } else if (i==1) {
 //       	for (int j=0; j<1000; j++) {
 //       		double correction = -22.6388 + 0.0246762*recoz[i].at(j) + 0.000818524*pow(recoz[i].at(j),2) - 1.32198*(1/pow(10,6))*pow(recoz[i].at(j),3) - 2.34171*(1/pow(10,9))*pow(recoz[i].at(j),4);
 //       		recoz[i].at(j) -= correction;
 //       		double all = -22.3876 + 0.0314674*recoz[i].at(j) + 0.000789842*pow(recoz[i].at(j),2) -1.47404*(1/pow(10,6))*pow(recoz[i].at(j),3) -1.84396*(1/pow(10,9))*pow(recoz[i].at(j),4);             
 //                       bz[i].at(j) = recoz[i].at(j) - all;
 //       	}
 //       } else if (i==2) {
 //       	for (int j=0; j<1000; j++) {
 //       		double correction = -23.184 - 0.00511548*recoz[i].at(j) + 0.000910545*pow(recoz[i].at(j),2) - 7.80652*(1/pow(10,7))*pow(recoz[i].at(j),3) - 3.88663*(1/pow(10,9))*pow(recoz[i].at(j),4);
 //       		recoz[i].at(j) -= correction;
 //       		double all = -22.9297 + 0.0229553*recoz[i].at(j) + 0.000804151*pow(recoz[i].at(j),2) -1.30832*(1/pow(10,6))*pow(recoz[i].at(j),3) -2.21049*(1/pow(10,9))*pow(recoz[i].at(j),4);             
 //                       bz[i].at(j) = recoz[i].at(j) - all;
 //       	}
 //       } else if (i==3) {
 //       	for (int j=0; j<1000; j++) {
 //       		double correction = -19.8847 - 0.00426127*recoz[i].at(j) + 0.000956135*pow(recoz[i].at(j),2) - 6.67826*(1/pow(10,7))*pow(recoz[i].at(j),3) - 4.70135*(1/pow(10,9))*pow(recoz[i].at(j),4);
 //       		recoz[i].at(j) -= correction;
 //       		double all = -19.355 + 0.020879*recoz[i].at(j) + 0.000809046*pow(recoz[i].at(j),2) -1.19208*(1/pow(10,6))*pow(recoz[i].at(j),3) -2.56328*(1/pow(10,9))*pow(recoz[i].at(j),4);             
 //                       bz[i].at(j) = recoz[i].at(j) - all;
 //       	}
 //       } else if (i==4) {
 //       	for (int j=0; j<1000; j++) {
 //       		double correction = -20.2096 - 0.00925059*recoz[i].at(j) + 0.000839815*pow(recoz[i].at(j),2) - 4.56888*(1/pow(10,7))*pow(recoz[i].at(j),3) - 3.77353*(1/pow(10,9))*pow(recoz[i].at(j),4);
 //       		recoz[i].at(j) -= correction;
 //       		double all = -18.5288 + 0.0419566*recoz[i].at(j) + 0.000614606*pow(recoz[i].at(j),2) -1.40293*(1/pow(10,6))*pow(recoz[i].at(j),3) -8.17669*(1/pow(10,10))*pow(recoz[i].at(j),4);             
 //                       bz[i].at(j) = recoz[i].at(j) - all;
 //       	}	
 //       } else if (i==5) {
 //       	for (int j=0; j<1000; j++) {
 //               	double correction = -25.1842 - 0.0127975*recoz[i].at(j) + 0.000952892*pow(recoz[i].at(j),2) - 2.5854*(1/pow(10,7))*pow(recoz[i].at(j),3) - 4.786258*(1/pow(10,9))*pow(recoz[i].at(j),4);
 //               	recoz[i].at(j) -= correction;
 //       		double all = -25.1892 + 0.0167561*recoz[i].at(j) +  0.000890128*pow(recoz[i].at(j),2) -7.96254*(1/pow(10,7))*pow(recoz[i].at(j),3) -3.653*(1/pow(10,9))*pow(recoz[i].at(j),4);             
 //                       bz[i].at(j) = recoz[i].at(j) - all;
 //       	}
 //       } else if (i==6) {
 //       	for (int j=0; j<1000; j++) {
 //               	double correction = -16.7996 + 0.0237262*recoz[i].at(j) + 0.000759375*pow(recoz[i].at(j),2) - 7.5112*(1/pow(10,7))*pow(recoz[i].at(j),3) - 3.27761*(1/pow(10,9))*pow(recoz[i].at(j),4);
 //       		recoz[i].at(j) -= correction;
 //       		double all = -18.6675 + 0.0379192*recoz[i].at(j) + 0.000781891*pow(recoz[i].at(j),2) -1.08859*(1/pow(10,6))*pow(recoz[i].at(j),3) -2.97615*(1/pow(10,9))*pow(recoz[i].at(j),4);             
 //                       bz[i].at(j) = recoz[i].at(j) - all;
 //       	}
 //       } else if (i==7) {
 //       	for (int j=0; j<1000; j++) {
 //               	double correction = -18.4332 - 0.0281257*recoz[i].at(j) + 0.000886945*pow(recoz[i].at(j),2) - 2.84511*(1/pow(10,7))*pow(recoz[i].at(j),3) - 4.24852*(1/pow(10,9))*pow(recoz[i].at(j),4);
 //       		recoz[i].at(j) -= correction;
 //               	double all = -19.1595 + -0.0169758*recoz[i].at(j) + 0.000881852*pow(recoz[i].at(j),2) -5.44949*(1/pow(10,7))*pow(recoz[i].at(j),3) -3.72381*(1/pow(10,9))*pow(recoz[i].at(j),4);             
 //                       bz[i].at(j) = recoz[i].at(j) - all;
 //       	}
 //       } else if (i==8) {
 //       	for (int j=0; j<1000; j++) {
 //               	double correction = -23.2351 - 0.00990602*recoz[i].at(j) + 0.000871257*pow(recoz[i].at(j),2) - 2.80623*(1/pow(10,7))*pow(recoz[i].at(j),3) - 4.02855*(1/pow(10,9))*pow(recoz[i].at(j),4);
 //       		recoz[i].at(j) -= correction;
 //       		double all = -22.8671 + 0.00179504*recoz[i].at(j) + 0.000814712*pow(recoz[i].at(j),2) -5.28401*(1/pow(10,7))*pow(recoz[i].at(j),3) -3.14937*(1/pow(10,9))*pow(recoz[i].at(j),4);             
 //                       bz[i].at(j) = recoz[i].at(j) - all;
 //       	}	
 //       } else if (i==9) {
 //       	for (int j=0; j<1000; j++) {
 //       		double correction = -17.4075 + 0.0113474*recoz[i].at(j) + 0.00079277*pow(recoz[i].at(j),2) - 5.148*(1/pow(10,7))*pow(recoz[i].at(j),3) - 3.22079*(1/pow(10,9))*pow(recoz[i].at(j),4);
 //       		recoz[i].at(j) -= correction;
 //               	double all = -16.6329 + 0.0203525*recoz[i].at(j) + 0.000723641*pow(recoz[i].at(j),2) -7.6036*(1/pow(10,7))*pow(recoz[i].at(j),3) -2.30981*(1/pow(10,9))*pow(recoz[i].at(j),4);             
 //                       bz[i].at(j) = recoz[i].at(j) - all;
 //       	}
 //       }
 // }


  //find mean reco error at each energy
  std::vector<double> zcorr_error[10];
  std::vector<double> mean_zcorr_error(10,0);
  for (int i=0; i<10; i++) {
        int size = true_vertex[i][0].size();
        for (int j=0; j<size; j++) {
                double xdiff = abs(reco_vertex[i][0].at(j)-true_vertex[i][0].at(j));
                double ydiff = abs(reco_vertex[i][1].at(j)-true_vertex[i][1].at(j));
                double zdiff = abs(reco_vertex[i][2].at(j)-true_vertex[i][2].at(j));
                zcorr_error[i].push_back(sqrt(pow(xdiff,2)+pow(ydiff,2)+pow(zdiff,2)));
                mean_zcorr_error[i] += zcorr_error[i].back();
        }
        mean_zcorr_error[i] = mean_zcorr_error[i]/size;
  }

  //find standard deviation at each energy
  std::vector<double> zcorr_sd(10,0);
  for (int i=0; i<10; i++) {
        int size = zcorr_error[i].size();
        for (int j=0; j<size; j++) {
                zcorr_sd[i] += pow(abs(zcorr_error[i].at(j)-mean_zcorr_error[i]),2);
        }
        zcorr_sd[i] = sqrt(zcorr_sd[i]/size);
  }



  //plot the mean new mean error at each energy level
  TGraphErrors *c3 = new TGraphErrors( 10, energies.data(), mean_zcorr_error.data(), blank.data(), zcorr_sd.data() );
  c3->SetMarkerColor(kGreen);
  c3->SetMarkerStyle(21);
  c3->SetName("corrected_reco_vertex_error");
  c3->GetXaxis()->SetTitle("Reco Z");
  c3->GetYaxis()->SetTitle("|Reco Z - True Z|");
  c3->Draw("apl");
  c3->Write();

//  mg->SetName("mean_mag_of_error");
//  mg->GetXaxis()->SetTitle("Energy MeV");
//  mg->GetYaxis()->SetTitle("Mean |Reco Z - True Z|");
//  mg->Draw("apl");
//  mg->Write();
//
  //gPad->Update();
  //gPad->Modified();

  //-----------------------------------------------------------------------------------------------------------------------------------------------------

  // make mean and sd variables
  //int energy_index = 0;
  //for (int i=0; i<10000;) {
  //      int nentries = 0;
  //      while ( nentries < 1000) {
  //      	mean_pmt_scint_hit_total.at(energy_index) += pmt_scint_hit_total.at(i);
  //      	mean_sipm_ceren_hit_total.at(energy_index) += sipm_ceren_hit_total.at(i);
  //      	nentries++;
  //      	i++;
  //      }
  //      mean_pmt_scint_hit_total.at(energy_index) = mean_pmt_scint_hit_total.at(energy_index)/nentries;
  //      mean_sipm_ceren_hit_total.at(energy_index) = mean_sipm_ceren_hit_total.at(energy_index)/nentries;
  //      energy_index++;
  //      //std::cout << i << std::endl;
  //}
 
 


 // OBstd::cout << "Made it out" << std::endl;

 // std::vector<double> energy_values = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50};

 // TGraphErrors *scint_yield = new TGraphErrors(10, energy_values.data(), mean_pmt_scint_hit_total.data());
 // scint_yield->SetMarkerColor(kRed);
 // scint_yield->SetMarkerStyle(21);
 // scint_yield->SetName("scint_yield");
 // scint_yield->Draw("apl");
 // scint_yield->Write();

 // TGraphErrors *ceren_yield = new TGraphErrors(10, energy_values.data(), mean_sipm_ceren_hit_total.data());
 // ceren_yield->SetMarkerColor(kBlue);
 // ceren_yield->SetMarkerStyle(21);
 // ceren_yield->SetName("ceren_yield");
 // ceren_yield->Draw("apl");
 // ceren_yield->Write();
 //



  // fill histograms
  //for (int i=0; i<mean_pmt_scint_hit_total.size(); i++) {
//	  std::cout << (5*i)+5 << " " << mean_pmt_scint_hit_total.at(i) << std::endl;
//	  hscintyield.Fill( (5*i)+5,  mean_pmt_scint_hit_total.at(i));
  //}
  //for (int i=0; i<mean_sipm_ceren_hit_total.size(); i++) {
//	  hcerenyield.Fill( (5*i)+5,  mean_sipm_ceren_hit_total.at(i));
  //}

  //hscintyield.Write();
  //hcerenyield.Write();

  fout->Close();
  
  return 0;
}
