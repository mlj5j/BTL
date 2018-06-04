#include <fstream>
#include <iostream>
#include <ios>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstring>
#include <string.h>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"


double findRMS(double ch[1024], float t[1024], float lowtime, float hightime);
double findmean(double ch[1024], float t[1024], float lowtime, float hightime);


//int convert_test(int argc, char **argv)
int convert_test(const char* name)
{

  //checking that there is a file specified for processing and then opening it.
  std::ifstream infile;
  infile.open(name, std::ios::binary | std::ios::in);
  std::cout << "Opening file " << name << " ...................." << std::endl;
  /*
    //This section was intended to make things easier and more general, but I couldn't get it to work...  
  if (argc == 2)
    {
      
      infile.open(argv[1], std::ios::binary | std::ios::in);
      //infile.open("SiPM-67V-Dark-RedAdvatech.dat", std::ios::binary | std::ios::in);
      std::cout << "Opening file " << argv[1] << " ............." << std::endl;
      std::cout << std::endl;
      if (!infile.is_open())
	{
	  std::cerr << "!! Error opening file: " << argv[1] << std::endl;
	  return 1;
	}
    }
  else
    {
      std::cerr << "!! There is no input file!!" << std::endl;
      return 1;
    }  
  int file_len = strlen(argv[1]); //determining length of the input file name 
  string filename = argv[1];  // initializing output filename to be same as input filename
  filename.replace(file_len -3, 3, "root"); //replacing "dat" at the end with "root"
*/

  //Creating a root file with the same name as the dat file
  int file_len = strlen(name);
  string filename = name;
  filename.replace(file_len -3, 3, "root");

  //Initializing the output ROOT file
  TFile *outfile = new TFile ((char *) filename.c_str (), "recreate");
  std::cout << "Creating ROOT file " << filename << " ......." << std::endl;
  std::cout << std::endl;

  //declaring variables to be used in ROOT file
  double ch1mV[1024]; //This is the voltage recorded on channel 1 in mV
  float time[1024];  //This is the time of a sample in ns with time[0] = 0
  int event = 0;
  double RMS;  //RMS of pedestal
  double amp;
  double mean;
  double integral;

  //Defining the function we will use to fit the pulses (Landau function)
  TF1 *fn1 = new TF1("ourfn", "landau(0)");

  //Initializing the Tree and branches in the ROOT file
  TTree *pulse = new TTree("pulse", "This is a tree");
  TBranch *ch1_b = pulse->Branch("ch1mV", &ch1mV, "ch1mV[1024]/D"); 
  TBranch *time_b = pulse->Branch("time", &time, "time[1024]/F");  
  TBranch *event_b = pulse->Branch("event", &event, "event/I");
  TBranch *RMS_b = pulse->Branch("RMS", &RMS, "RMS/D");
  TBranch *amp_b = pulse->Branch("amp", &amp, "amp/D");
  TBranch *mean_b = pulse->Branch("mean", &mean, "mean/D");
  TBranch *integral_b = pulse->Branch("integral", &integral, "integral/D");

  //Declaring some dummy variables for reading headers we don't care about.
  char tmpHeader[4];
  char tmpsmallHeader[2];
  int SerialNumber;
  short Date[8];
  unsigned short ch1[1024]; //This will be the raw value of the channel, which we will convert to mV later.
  bool endoffile = false;
  char tmpBoardNumber[5];

  //Reading bytes until we get to stuff we actually care about.  
  infile.read((char *) &tmpHeader, 4);
  std::cout << "header 1 reads " << tmpHeader << std::endl;
  infile.read((char *) &tmpHeader, 4);
  std::cout << "header 2 reads " << tmpHeader << std::endl; 
  infile.read((char *) &tmpsmallHeader, 2);
  infile.read((char *) &SerialNumber, 2);
  std::cout << "header 3 reads " << tmpsmallHeader << " " << SerialNumber << std::endl;
  infile.read((char *) &tmpHeader, 4);
  std::cout << "header 4 reads " << tmpHeader << std::endl;

  infile.read((char *) &time, 4096);  //This gets the time values in nanoseconds for each sample in the buffer.
  for (int i=1; i<1024; i++)
    {
      time[i] += time[i-1];
    }
  std::cout << "Time calculated" << std::endl;
  
  while(!infile.eof())
    {
      if (event%500==0)
	std::cout << "Processing event #" << event << std::endl;

      //Reading more bytes that we don't care about...
      infile.read((char *) &tmpHeader, 4);
      //std::cout << "beginning of loop: " << tmpHeader << std::endl;
      infile.read((char *) &SerialNumber, 4);
      infile.read((char *) &Date, 16);
      infile.read((char *) &tmpBoardNumber, 4);
      infile.read((char *) &tmpBoardNumber, 4);
      infile.read((char *) &tmpHeader, 4);
      infile.read((char *) &SerialNumber, 4);
      // std::cout << "Scaler # " << SerialNumber << std::endl;

      //Reading the values of each sample in the buffer for an event
      infile.read((char *) &ch1, 2048);

      //Converting raw value to mV
      for (int j=0; j<1024; j++)
	{
	  ch1mV[j] = ch1[j]/65535.-0.5;

	}

      TGraph *g1 = new TGraph();  //Making a TGraph for each event that we will use to fit each pulse with a function

      RMS = findRMS(ch1mV, time, 40, 120);
      mean = findmean(ch1mV, time, 40, 120);

      //Filling TGraph for each event (pulse)
      for (int j=0; j<1024; j++)
	{
	  ch1mV[j] = mean - ch1mV[j];
	  g1->SetPoint(j, time[j], ch1mV[j]);
	  
	}
     
      amp = -ch1mV[TMath::LocMin(1024, ch1mV)];  //Finding amplitude of the pulse

      fn1->SetRange(50, 110);  //Setting fit range the match the region in which we see most pulses
      
      g1->Fit(fn1, "QR"); //Fitting without printing out all of the statistics to the screen (Quiet mode)
      
      //Drawing a pulse to check that the fit is working
      if (event==100)
	{
	  TCanvas *c2 =  new TCanvas("c2", "c2", 500, 500);
	  g1->Draw();
	}
      //Integrating the fit function to get the area under the pulse (related to net charge of the pulse)
      TF1 *fn1 = g1->GetFunction("ourfn");
      integral = fn1->Integral(50, 110);

      event++;
      pulse->Fill(); //Filling the tree's branches
      if (infile.eof())
        {
	  std::cout << "End of file reached ......." << std::endl;
	  std::cout << "Processed " << event << " events ..." << std::endl;
          endoffile = true;
          break;
        }

    }
  std::cout << "Tree was filled" << std::endl;
  outfile->Write();
  outfile->Close();
  std::cout << "ROOT file written and closed" << std::endl;

  infile.close();

  return 0;

}
double findmean(double ch[1024], float t[1024], float lowtime, float hightime)
{
  double sum = 0;
  int counter = 1;
  for (int i=0; i<1024; i++)
    {
      if (t[i] < lowtime || t[i] > hightime)
	{
	  sum += ch[i];
	  counter++;
	  
	}
      
    }
  double mean = sum/counter;
  return mean;
}

double findRMS(double ch[1024], float t[1024], float lowtime, float hightime)
{
  double sum = 0;
  int counter = 1;
  for (int i=0; i<1024; i++)
    {
      if (t[i] < lowtime || t[i] > hightime)
	{
	  sum += ch[i]*ch[i];
	  counter++;
	  
	}
      
    }
  double RMS = TMath::Sqrt(sum/counter);
  return RMS;
}
