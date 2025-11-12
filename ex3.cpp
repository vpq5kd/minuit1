// A simple C++ program to illustrate the use of ROOT class TMinuit 
// for function minimization.  The example shows a Maximum Likelihood
// fit for the mean of an exponential pdf in which TMinuit 
// minimizes - 2*log(L).   The user must specify what to minimize in the 
// function fcn, shown in this file below.

// fcn passes back f = -2*ln L by reference; this is the function to minimize.
// The factor of -2 allows MINUIT to get the errors using the same
// recipe as for least squares, i.e., go up from the minimum by 1.

// TMinuit does not allow fcn to be a member function, and the function
// arguments are fixed, so one of the only ways to bring the data  
// into fcn is to declare a pointer to the data (xVecPtr) as global.

// For more info on TMinuit see root.cern.ch/root/html/TMinuit.html .

// Based on example by 
// Glen Cowan
// RHUL Physics
// 4 December 2006

// Update Bob Hirosky: Sep 2013


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include <TMinuit.h>
#include <TRandom2.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>
#include <TAxis.h>
#include <TLine.h>
#include <TFile.h>
#include <TLatex.h>
#include <TH2F.h>
#include <TF2.h>
using namespace std;

// ---------- GLOBALS ----------
// Declare pointer to data as global
// The use of global variables is a disfavored coding practice, but
// in this case it will allow us to avoid writing more complex code
// also passing the fit function, to make the objective fcn more generic
TF2 *fparam;

 TFile *file = TFile::Open("fitInputs.root", "READ");
  TH2F  *hdata = (TH2F*)file->Get("hdata");
  TH2F  *hbkg = (TH2F*)file->Get("hbkg");
  
//-------------------------------------------------------------------------
// The pdf to be fitted, here an exponential.
// This is a convenient interface for associating the model
// with a TF1 later
//-------------------------------------------------------------------------
double gauss2D(double *xyPtr, double par[]){
	double x = xyPtr[0];
	double y = xyPtr[1];
	
	double A = par[0];
	double mux = par[1];
	double sigmax = par[2];
	double muy = par[3];
	double sigmay = par[4];
	double cbkg = par[5];

	double signal = A * exp(-pow((x-mux)/sigmax, 2))
			  * exp(-pow((y-muy)/sigmay, 2));
	double bkg = cbkg *hbkg->Interpolate(x,y);
	return signal + bkg;
}
// This is our OBJECTIVE Function
// return NLL given a histogram and function

double calcNLL(TH1F* h, TF1* f){
  double nll=0;
  for (int i=1; i<=h->GetNbinsX(); i++){
    double x=h->GetBinCenter(i);
    int n=(int)(h->GetBinContent(i));
    double mu=f->Eval(x);
    if (mu<1e-10) mu=1e-10;    // avoid log(0) problems if we go outside a reasonable range!
    nll -= n * TMath::Log(mu) - mu  - TMath::LnGamma(n+1);
  }
  // cout << "nll "<< nll <<endl;
  return 2*nll;   // factor of 2 so the 1 sigma error contours follow the chi^2 convention
}

//utilized chatGPT to aid with creating the chi^two function
double calcChi2(TH2F *h, TF2 *f){
	double chi2 = 0.0;
	for (int ix = 1; ix <= h->GetNbinsX(); ix++){
		for (int iy = 1; iy <= h->GetNbinsY(); iy++){
			double x = h->GetXaxis()->GetBinCenter(ix);
			double y = h->GetYaxis()->GetBinCenter(iy);

			double obs = h->GetBinContent(ix, iy);
			double err = h->GetBinError(ix, iy);
			if (err <= 0) err = sqrt(obs +1);

			double expval = f->Eval(x,y);
			chi2 += pow((obs - expval)/err,2);	
		}
	}
	return chi2;
}

//-------------------------------------------------------------------------

void fcn(int& npar, double* deriv, double& f, double par[], int flag){

  for (int i=0; i<npar; i++){
    fparam->SetParameter(i,par[i]);
  }

  f = calcChi2(hdata,fparam);
 
}

int main(int argc, char **argv) {

  // This allows you to view ROOT-based graphics in your C++ program
  // If you don't want view graphics (eg just read/process/write data files), 
  // this line can be removed
  TApplication theApp("App", &argc, argv);

  TCanvas* canvas = new TCanvas();
  // ***************************************
  // not important for this example
  // Set a bunch of parameters to make the plot look nice
  canvas->SetFillColor(0);
  canvas->UseCurrentStyle();
  canvas->SetBorderMode(0);        
  canvas->SetFrameBorderMode(0);  
  //gROOT->SetStyle("Plain");
  canvas->UseCurrentStyle();
  gROOT->ForceStyle();

  //gStyle->SetOptStat(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleSize(0.04);
  gStyle->SetTitleFont(42, "hxy");      // for histogram and axis titles
  gStyle->SetLabelFont(42, "xyz");      // for axis labels (values)
  gROOT->ForceStyle();
  // ***************************************
 const int npar = 6;
  double xmin = hdata->GetXaxis()->GetXmin();
  double xmax = hdata->GetXaxis()->GetXmax();
  double ymin = hdata->GetYaxis()->GetXmin();
  double ymax = hdata->GetYaxis()->GetXmax();
  TF2* fitfunc  = new TF2("fitfunc", gauss2D, xmin, xmax, ymin, ymax, npar);
  TString names[npar] = {"A", "mux", "sigx", "muy","sigy", "cbkg"};
  double par[npar] = {50, 3.4, 1.2, 2.0, 1.4, .4};
  double stepSize[npar] = {1,.1,.1,.1,.1,.01};

  TMinuit minuit1(npar);
  minuit1.SetFCN(fcn);

  fparam=fitfunc;  // our model

  for (int i = 0; i<npar; i++){
  	minuit1.DefineParameter(i, names[i].Data(), par[i], stepSize[i], 0, 0);
  }
  // Do the minimization!
  minuit1.Migrad();       // Minuit's best minimization algorithm

  // Get the result
  double outpar1[npar], err1[npar];
  for (int i=0; i<npar; i++){
    minuit1.GetParameter(i,outpar1[i],err1[i]);
  };
  fparam->SetParameters(outpar1);

  double fmin2, fedm2, errdef2; 
  int npari2, nparx2, istat2;
  minuit1.mnstat(fmin2, fedm2, errdef2, npari2, nparx2, istat2);

  TCanvas *tc = new TCanvas("tc", "Canvas", 1200, 900);
  tc->Divide(2,2);

  tc->cd(1);
  hdata->Draw("lego2");

  tc->cd(2);
  TH2F *hfit = (TH2F*)hdata->Clone("hfit");
  hfit->Reset();
  for(int ix = 1; ix<= hfit->GetNbinsX(); ix++){
  	for (int iy = 1; iy <= hfit->GetNbinsY(); iy++){
		double x = hfit->GetXaxis()->GetBinCenter(ix);
		double y = hfit->GetYaxis()->GetBinCenter(iy);
		double val = fparam->Eval(x,y);
		hfit->SetBinContent(ix,iy,val);
	}
  }
  hfit->SetTitle("Fit Results (Signal + Background)");
  hfit->Draw("lego2");

  tc->cd(3);
  TH2F *hres = (TH2F*)hdata->Clone("hres");
  for (int ix = 1; ix <= hres->GetNbinsX(); ix++){
  	for (int iy = 1; iy<=hres->GetNbinsY(); iy ++){
		double difference = hdata->GetBinContent(ix, iy) - hfit->GetBinContent(ix,iy);
		hres->SetBinContent(ix, iy, difference);
	}
  
  }
  hres->SetTitle("Residuals");
  hres->Draw("lego2");


  tc->cd(4);
  TH2F *hbkg_best = (TH2F*)hbkg->Clone("hbkg_best");
  hbkg_best->Scale(fparam->GetParameter(5));
  TH2F *hsub = (TH2F*)hdata->Clone("hsub");
  hsub->Add(hbkg_best, -1);
  hsub->SetTitle("Data After Subtracting the Best Fit Background");
  hsub->Draw("lego2");

  tc->SaveAs("ex3_fig.png");

  for (int i = 0; i<npar;i++){
  	std::cout << names[i] << "=" << outpar1[i] << "+/-" << err1[i] << std::endl;
  }
  theApp.Run(true);
  canvas->Close();
  
  return 0;

}
