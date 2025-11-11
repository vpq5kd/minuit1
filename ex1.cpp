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
using namespace std;

// ---------- GLOBALS ----------
// Declare pointer to data as global
// The use of global variables is a disfavored coding practice, but
// in this case it will allow us to avoid writing more complex code
TH1F *hdata;
// also passing the fit function, to make the objective fcn more generic
TF1 *fparam;

//-------------------------------------------------------------------------
// The pdf to be fitted, here an exponential.
// This is a convenient interface for associating the model
// with a TF1 later
double doubleGauss(double *xPtr, double par[]){
	double x = *xPtr;
	double A1 = par[0];
	double mu1 = par[1];
	double sigma1 = par[2];
	double A2 = par[3];
	double mu2 = par[4];
	double sigma2 = par[5];
	return A1 * TMath::Gaus(x,mu1, sigma1, true)+
		A2 * TMath::Gaus(x, mu2, sigma2, true);
}

double gumbelPDF(double *xPtr, double par[]){
	double x = *xPtr;
	double A = par[0];
	double mu = par[1];
	double beta = par[2];
	double z = (x-mu)/beta;
	return A * (1.0/beta) *exp(-(z + exp(-z)));

}

//-------------------------------------------------------------------------
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
double calcChi2(TH1F *h, TF1 *f){
	double chi2 = 0.0;
	for (int i = 1; i <= h->GetNbinsX(); i++){
		double x = h->GetBinCenter(i);
		double y_o = h->GetBinContent(i);
		double y_e = f->Eval(x);
		double err = h->GetBinError(i);

		chi2+= pow((y_o - y_e)/err,2);
	}
	return chi2;
}

//-------------------------------------------------------------------------
// Minuit fcn: calculates value of the function to be minimized using
// the data and the model function
// This is the interface used to define our objective function.  We use
// only a subset of the input parameters below.
// npar: number of parameters
// par: array of parameter values
// f: the value of the objective function
// Minuit can also pass the gradient(deriv) of the objective function wrt the
// current parameters or flags that can trigger special operations (eg perform some initialization)

void fcn(int& npar, double* deriv, double& f, double par[], int flag){

  for (int i=0; i<npar; i++){
    fparam->SetParameter(i,par[i]);
  }

  f = calcChi2(hdata,fparam);
 
}

//-------------------------------------------------------------------------

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
  TFile *file = TFile::Open("distros.root", "READ");
  TH1F  *hexp = (TH1F*)file->Get("dist1");
  double xmin = hexp->GetXaxis()->GetXmin();
  double xmax = hexp->GetXaxis()->GetXmax();

  
  const int npar1 = 6;
  TF1* gauss2 = new TF1("gauss2", doubleGauss, xmin, xmax, npar1);
  //used chatGPT to help with parameter guessing
  double par1[npar1] = {
  	hdata->GetMaximum(), hdata->GetMean() -1, 0.5,
	hdata->GetMaximum()/2, hdata->GetMean() +1, 0.5
  };
  double stepSize1[npar1] = {.1,.1,.1,.1,.1,.1};
  double minVal1[npar1] = {0,0,0,0,0,0};
  double maxVal1[npar1] = {0,0,0,0,0,0};
  TString parName1[npar1] = {"A1", "mu1", "sigma1", "A2","mu2","sigma2"};

  TMinuit minuit1(npar1);
  minuit1.SetFCN(fcn);

  hdata=hexp;     // histogram to fit
  fparam=gauss2;  // our model

  for (int i = 0; i<npar1; i++){
  	minuit1.DefineParameter(i, parName1[i].Data(), par1[i], stepSize1[i], minVal1[i], maxVal1[i]);
  }
  // Do the minimization!
  minuit1.Migrad();       // Minuit's best minimization algorithm

  // Get the result
  double outpar1[npar1], err1[npar1];
  for (int i=0; i<npar1; i++){
    minuit1.GetParameter(i,outpar1[i],err1[i]);
  };
  gauss2->SetParameters(outpar1);

  gauss2->SetLineStyle(1);             //  1 = solid, 2 = dashed, 3 = dotted
  gauss2->SetLineColor(1);             //  black (default)
  gauss2->SetLineWidth(2);

  gauss2->GetXaxis()->SetTitle("x");
  gauss2->GetYaxis()->SetTitle("f(x; fit from #chi^{2})");

  // Plot the result.
  hdata->Draw("e");
  gauss2->Draw("same");

  // summarize the fitting results
  cout << "\n==========================\n"<<endl;
  double fmin, fedm, errdef;
  int npari, nparx, istat;  // see https://root.cern/doc/master/classTMinuit.html 
  minuit1.mnstat(fmin, fedm, errdef, npari, nparx, istat);
  cout << "minimum of chi^2 = " << fmin << endl;
  cout << "fit status = " << istat << endl;
  cout << "best fit parameters\n" <<endl;
  for (int i=0; i<npar1; ++i){
    cout << i << " : " << outpar1[i] << " +- " << err1[i] << endl;
  };

  cout << "\nTo exit, quit ROOT from the File menu of the plot (or use control-C)" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program
  theApp.Run(true);
  canvas->Close();

  return 0;

}
