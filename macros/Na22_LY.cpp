#include <iostream>
#include <string.h>

#include "TProfile.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TSpectrum2.h"
#include "TH2.h"
#include "TH1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TBranch.h"

int FindLeadingEdge(double ch[1024], float t[1024], double th);

double GetBaseline(double ch[1024], float t[1024], int startbin, int nsamples);

double GetStdDev(double ch[1024], double mean);

void Na22_LY(const char* name)
{
  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1);

  double chmV[1][1024];
  float time[1][1024];
  double amp[1];
  double threshold = 0.1; //This is threshold voltage in mV
  int tstartbin = 0;
  double baseline = -50;
  Int_t nfound = 0;
  float tot = 0;
  
  string filename = name;
  TFile *f = new TFile((char *) filename.c_str ());
  TTree *pulse = (TTree*)f->Get("pulse");

  pulse->SetBranchAddress("chmV", &chmV);
  pulse->SetBranchAddress("time", &time);
  pulse->SetBranchAddress("amp", &amp);
  TH1F *hsigma = new TH1F("hsigma", "Noise (stddev)", 50, 0, 4);
  TH1F *hbaseline = new TH1F("hbaseline", "Baseline", 200, -6, 16);
  TH1F *hamp = new TH1F("hamp", "Pulse Amplitude", 50, 0, 10);
  TH1F *hQ = new TH1F("hQ", "Integrated Charge", 100,0, 200);
  TH1F *hrt = new TH1F("hrt", "Risetime", 100, 0, 18);
  
  TProfile *peakcheck = new TProfile("peakcheck", "Peak profiles", 1024, 0, 200);
  TGraph *g1 = new TGraph();
  TGraph *g2 = new TGraph();
  TGraph *g3 = new TGraph();
  
  bool foundlow, foundhigh;
  int lowbin, highbin;

  int file_len = strlen(name);
  string fileoutname = name;
  fileoutname.insert(file_len -5, "_output");
  TFile *fout = new TFile((char *) fileoutname.c_str (), "RECREATE");
  if (fout->IsOpen())
    std::cout << "Output file ready for data" << std::endl;
  TTree *tree = new TTree("tree", "LY data tree");

  double sigma;
  double charge;
  float risetime;
  double chmax;

  //Defining branches of output file
  TBranch *ch_b = tree->Branch("chmV", &chmV, "chmV[1][1024]/D");
  TBranch *time_b = tree->Branch("time", &time, "time[1][1024]/F");

  TBranch *risetime_b = tree->Branch("risetime", &risetime, "risetime/F");
  TBranch *sigma_b = tree->Branch("sigma", &sigma, "sigma/D");
  TBranch *charge_b = tree->Branch("charge", &charge, "charge/D");
  TBranch *chmax_b = tree->Branch("chmax", &chmax, "chmax/D");
  TBranch *tot_b = tree->Branch("tot", &tot, "tot/F");
  
  for (int i=0; i<pulse->GetEntries(); i++)
    {
      //Initialize values to be recorded
      sigma = 0;
      charge = 0;
      risetime = 0;
      chmax = 0;
      //Initialize flags stating if the beginning and ending bins of the pulse have been found
      foundlow = false;
      foundhigh = false;
      
      if(i%1000==0)
	{
	  std::cout<< "processing event " << i << std::endl;
	}
      pulse->GetEntry(i);
      baseline = GetBaseline(chmV[0], time[0], 25, 60);
      sigma = GetStdDev(chmV[0], baseline);

      //for (int n=0; n<1024; n++)
      //chmV[0][n] -= baseline;
      
      hbaseline->Fill(sigma);
      tstartbin = FindLeadingEdge(chmV[0], time[0], 0.1 + baseline);
      if (tstartbin<50)
	continue;
      
      //  tstartbin -= 20;
      
      
      TProfile *ptemp = new TProfile("ptemp", "Temporary Profile", 1024, 0, 200);
      //Loope through the buffer starting with the first bin of the pulse
      for (int j=tstartbin; j<1024; j++)
	{
	  //If we haven't tagged a lowbin and the sample is above threshold, then tag this one as the lowbin
	  if (!foundlow && chmV[0][j] > (threshold - baseline))  //Changed!!!!!!!!!!!!!
	    {
	      lowbin = j;
	      foundlow = true;
	    }
	  //If we have tagged a lowbin, haven't found the highbin, and the sample is below the negative of the threshold, the tage this one as the highbin.
	  //Using the negative here is just a way of looking for a switch in polarity.  Maybe changing "-threshold" to zero would make more sense?
	  if (!foundhigh && foundlow && chmV[0][j] < -(threshold-baseline)) //Changed!!!!!!!!!!!!!!!!
	    {
	      highbin = j;
	      foundhigh = true;
	    }
	  ptemp->Fill(time[0][j], chmV[0][j]);
	  //	  peakcheck->Fill(time[0][j-tstartbin], chmV[0][j]);
	}
      //If the pulse is wide enough but not too wide, integrate it.  
      if (highbin-lowbin>80 && highbin-lowbin<800)
	{
	  charge = ptemp->Integral(lowbin, highbin-1)*(time[0][highbin]-time[0][lowbin])/(highbin-lowbin);
	  
	  for (int r=lowbin; r<highbin; r++)
	    {
	      if(chmV[0][r]>chmax)
		{
		  chmax=chmV[0][r];
		  risetime = time[0][r]-time[0][lowbin];
		}
	      peakcheck->Fill(time[0][r-lowbin], chmV[0][r]);
	    }
	  
	  g3->SetPoint(nfound,sigma, charge);  
	  nfound++;
	  if (nfound==1)
	    {
	      //cout << "lowbin is " << lowbin << endl;
	      // cout << "highbin is "<< highbin << endl;
	      for (int w=0; w<1024; w++)
		{
		  g1->SetPoint(w, time[0][w], chmV[0][w]);
		}
	      for (int w=lowbin; w<highbin; w++)
		{
		  g2->SetPoint(w-lowbin, time[0][w], chmV[0][w]);
		}
	   
	    }
	  tot = time[0][highbin] - time[0][lowbin];
	  hamp->Fill(chmax);
	  hrt->Fill(risetime);
	  hQ->Fill(ptemp->Integral(lowbin, highbin-1)*(time[0][highbin]-time[0][lowbin])/(highbin-lowbin));
	  tree->Fill();
	  if(nfound==1)
	    cout << "Charge is " << ptemp->Integral(lowbin, highbin-1)*(time[0][highbin]-time[0][lowbin])/(highbin-lowbin) << endl;
	  if (ptemp->Integral(lowbin, highbin-1)*(time[0][highbin]-time[0][lowbin])/(highbin-lowbin) > 130 && ptemp->Integral(lowbin, highbin-1)*(time[0][highbin]-time[0][lowbin])/(highbin-lowbin) < 150)
	    cout << "event number is " << i << endl;
	}
      delete ptemp;
    }

  TCanvas *cQ = new TCanvas("cQ", "Charge", 500, 500);
  hQ->Draw();

  TCanvas *cb = new TCanvas("cb", "Baseline", 500, 500);
  hbaseline->Draw();

  g2->SetLineColor(kRed);
  TCanvas *cprof = new TCanvas("cprof", "Profile", 500, 500);
  g1->Draw("Al");
  g2->Draw("lsame");

  TCanvas *c3 = new TCanvas("c3", "Sigma vs charge", 500, 500);
  g3->Draw("Ap");

  TCanvas *camp = new TCanvas("camp", "Amplitude", 500, 500);
  hamp->Draw();

  TCanvas *crt = new TCanvas("crt", "Risetime", 500, 500);
  hrt->Draw();

  fout->Write("tree");
  fout->Close();
  // f->Close();
  
}

//This function loops through the samples in the buffer until it finds a pulse that's at least 10ns wide.
int FindLeadingEdge(double ch[1024], float t[1024], double th)
{ 
  bool trig1 = false;
  float t1 = 0;
  int t1bin = 0;
  float t2 = 1;
  float ToT = 0;
  //loop through each sample in the buffer
  for(int k=0; k<1024; k++)
    {
      //Find the first sample in the buffer above the threshold voltage
      if (ch[k]>th && (!trig1))
	{
	  t1 = t[k]; //Time in ns of first sample in pulse
	  trig1=true;
       	  t1bin = k; //Bin number (or sample number) for first sample in pulse
	}
      //Once the first sample above threshold is found, find the next sample to go below threshold
      if (ch[k]<th && trig1)
	{
	  t2 = t[k];  //Time in ns of first sample after the pulse
	  ToT = t2 - t1;  //Width of the pulse in ns
	  trig1=false; //Set lowbin trigger to be false so if the loops isn't broken, we will look for another lowbin
	  //If the pulse width is more than 10ns we will break the loop which will stop the search for a pulse 
	  if (ToT > 10)
	    {
	      //	      t1bin = k;  //Bin number (or sample number for first sample in pulse
	      break;
	    }
	}
    }
  
  
  return t1bin;
}


double GetBaseline(double ch[1024], float t[1024], int startbin, int nsamples)
{
  double bline = -100;
  TProfile *pb = new TProfile("pb", "pb", 1024, 0, 200);
  for (int k=0; k<1024; k++)
    {
      pb->Fill(t[k], ch[k]);
    }
  bline = pb->Integral(startbin, startbin + nsamples -1)/nsamples;
  delete pb;
  return bline;
  
}

double GetStdDev(double ch[1024], double mean)
{
  double stdDev = 0;
  double sum = 0;
  for (int k=0; k<1024; k++)
    {
      sum += TMath::Power(ch[k] - mean,2);
    }
  stdDev = TMath::Sqrt(sum/1024);
  return stdDev;
}
