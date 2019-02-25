#include <iostream>

#include "TSpectrum.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"

void SinglePEfinder(const char* name)
{
  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1);

  //Declare variables
  double charge;

  //Load file and Tree
  string filename = name;
  TFile *f = new TFile((char *) filename.c_str());
  TTree *tree = (TTree*)f->Get("tree");
  

  //Set branch addresses
  tree->SetBranchAddress("charge",&charge);

  //Make histogram(s)
  TH1F *hQ = new TH1F("hQ", "Charge", 100, 0, 300);
  TGraph *gPeaks = new TGraph();
  //Loop through events
  for (int i=0; i<tree->GetEntries(); i++)
    {
      if (i%1000==0)
	std::cout<< "Processing Event " << i << std::endl;
      
      tree->GetEntry(i);
      hQ->Fill(charge);
      
    }
  TSpectrum *peakFinder = new TSpectrum(6);
  int Npeaks = peakFinder->Search(hQ, 2, "", 0.05);

  double* peak;
  peak = new double[Npeaks];
  peak = peakFinder->GetPositionX();

  TCanvas *cQ = new TCanvas("cQ", "Charge", 500, 500);
  gPad->SetLogy();
  hQ->Draw();
  for (int j=0; j<Npeaks; j++)
    {
      gPeaks->SetPoint(j,(double)j + 1.0, peak[j]);
    }
  TCanvas *cPeaks = new TCanvas("cPeaks", "Pe Peaks", 500, 500);
  gPeaks->Draw("A*");
 
}
