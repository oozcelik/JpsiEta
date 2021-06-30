#include <TROOT.h>
#include "TMath.h"
#include <iostream>

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "TPaveText.h"
#include "TLegend.h"
#include <TFile.h>

#include <TTree.h>

#include "RooNumIntConfig.h"
#include <RooRandom.h>

#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLatex.h"
#include "RooHist.h"

#include "RooExponential.h"
#include "RooGaussModel.h"
#include "TCanvas.h"
#include "RooLognormal.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TCut.h"
#include "RooGenericPdf.h"
#include "RooTruthModel.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooProdPdf.h"
#include "RooEffProd.h"
#include <TRandom3.h>
#include "RooRandom.h"
#include "RooPoisson.h"
#include "RooGamma.h"
#include "RooBernstein.h"
#include "RooGaussian.h"
#include "iostream"
#include "fstream"
#include "RooFitResult.h"
#include "RooCBShape.h"
#include "RooConstVar.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"


void trecoe(std::string filename="h.root"){

  TFile *theFile = new TFile(filename.c_str());
  TTree *tree = (TTree *)theFile->Get("DecayTree");

  using namespace std;
  using namespace RooFit;  
//  RooRandom::randomGenerator()->SetSeed(1234);
  TRandom3 ran(0);
  TH1D *h1 = new TH1D("h1", "h1", 100, -0.5, 0.5);
 
  RooRealVar* time =new RooRealVar("time","t_{reco}[ps]", -0.3, 0.3) ;
  //RooRealVar tresoe("tresoe", "tresoe", -10., 10.); 

  double dte;     tree->SetBranchAddress("B_LOKI_DTF_CTAUERR", &dte);
  double tresoe;

  int nevent = tree->GetEntries();

   for (int i = 0; i < nevent; ++i)
    {
         tree->GetEntry(i);
         tresoe = dte/0.299792458;
     for (int j=0; j<5; j++){

         double var  = ran.Gaus(0, tresoe);
         h1->Fill(var);
      //   cout << ran.GetSeed() << endl;
      }
   }
//  h1->Draw();

   RooDataHist datahist("hist", "hist", *time, h1) ;
/*   RooPlot *xframe = time->frame(100);
   TCanvas * c=new TCanvas("c","c",800,600);
   datahist.plotOn(xframe);
    xframe->Draw() ;*/


  // resolution model //

  RooRealVar gmm ("gmm", "gmm", 0);
  RooRealVar* gms1=new RooRealVar("gms1","gms1",0.03, -1., 1.);
  RooRealVar* gms2=new RooRealVar("gms2","gms2",0.1, 0, 1.2);
  RooRealVar* gms3=new RooRealVar("gms3","gms3",0.05, 0, 1.);


  RooRealVar* trecoe=new RooRealVar("trecoe","terr[ps]",0.003,0.26); // decay time per-event error. 
  RooTruthModel* CM =new RooTruthModel("CM","",*time); // to be added!
  
  RooGaussian *gausspdf1 = new RooGaussian("gausspdf1", "gausspdf1",  *time, gmm , *gms1) ;
  RooGaussian *gausspdf2 = new RooGaussian("gausspdf2", "gausspdf2",  *time, gmm , *gms2 ) ;
  RooGaussian *gausspdf3 = new RooGaussian("gausspdf3", "gausspdf3",  *time, gmm , *gms3 ) ;

   RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.6,0.,1.) ;
   RooRealVar sig2frac("sig2frac","fraction of component 2 in signal",0.8,0.,1.) ;
   RooAddPdf gausspdf("sig","Signal",RooArgList(*gausspdf1, *gausspdf2, *gausspdf3),RooArgList(sig1frac, sig2frac)) ;

    gausspdf.fitTo(datahist, RooFit::Optimize(false));
    RooPlot *xframe = time->frame(100);
    datahist.plotOn(xframe);
//    model.plotOn(xframe, Components(gausspdf), LineStyle(ELineStyle::kDashed), LineColor(kCyan));
  //  model.plotOn(xframe, Components(expoPdf),  LineStyle(ELineStyle::kDashed), LineColor(kGreen));
//    model.plotOn(xframe, Components(wpvPdf), LineStyle(ELineStyle::kDashed), LineColor(kRed));
    gausspdf.paramOn(xframe, Format("NE",FixedPrecision(5)));
    gausspdf.plotOn(xframe);
  
    TCanvas * c22 = new TCanvas();
//    c22->SetLogy();
    xframe->Draw();
    xframe->SetTitle("");

    gPad->SetBottomMargin(0.36);

    RooHist *xpull = xframe->pullHist();
    RooPlot *Pullframe = time->frame();

    Pullframe->addPlotable(xpull, "P");

    Pullframe->GetXaxis()->SetTitle("trecoE");
    Pullframe->GetYaxis()->SetNdivisions(205, kTRUE);
    Pullframe->GetYaxis()->SetTickLength(0.083);
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 1.0 - 0.06, 0, 0, 0.);
    pad2->SetTopMargin(1.0 - 0.36);
    pad2->SetFillStyle(0);
    pad2->Draw();
    pad2->cd();
    Pullframe->Draw();
}

