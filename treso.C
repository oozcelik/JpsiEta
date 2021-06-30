#include <TROOT.h>
#include "TMath.h"
#include <iostream>

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"
#include "RooFFTConvPdf.h"

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
#include "RooAbsCachedPdf.h"
#include "RooRealProxy.h"
#include "RooSetProxy.h"
#include "RooAbsReal.h"
#include "RooHistPdf.h"
#include "/Users/ozlemozcelik/work/txtfiles/TVirtualFFT.h"

void treco(std::string filename="h.root"){

  TFile *theFile = new TFile(filename.c_str());
  TTree *tree = (TTree *)theFile->Get("DecayTree");

  using namespace std;
  using namespace RooFit;  

  RooRealVar* time =new RooRealVar("time","t_{reco}[ps]", -0.3, 1.5) ;

  RooRealVar tau0( "tau0", "tau0", 1.2, -0.2, 2.);
  RooRealVar tau( "tau", "tau", 1.58, -0.2, 2.);
  RooTruthModel tm("tm","",*time); // Delta function

  RooRealVar scale ("scale", "scale", 0., 10.);

  // Gaussian resolution model //
  RooRealVar* gms1=new RooRealVar("gms1", "gms1", 0.031, 0., 0.05);
  RooRealVar* gms2=new RooRealVar("gms2","gms2", 0.05390, 0., 0.06);
  RooRealVar* gms3=new RooRealVar("gms3","gms3",0.0898, 0.06, 0.1);

  RooRealVar gmm ("gmm", "gmm", -0.004, -0.2, 1.);
  RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.362) ;
  RooRealVar sig2frac("sig2frac","fraction of component 2 in signal",0.523) ;

  RooGaussModel gm1( "gm1", "gauss model-1", *time, gmm, *gms1);
  RooGaussModel gm2( "gm2", "gauss model-2", *time, gmm, *gms2);
  RooGaussModel gm3( "gm3", "gauss model-3", *time, gmm, *gms3);

   RooAddModel gaussmodel( "gaussmodel","gaussmodel",RooArgList(gm1, gm2, gm3), RooArgList(sig1frac, sig2frac));

 //  RooFFTConvPdf conv("conv","delta (X) gaussian", *time, tm, gaussmodel); 

  /// effective resolution 
   RooRealVar effs1 ("effs1", "effs1", 0.052864652, 0., 0.5);
   RooGaussModel effgaussmodel( "effgm", "effgm", *time, gmm, effs1);

  //Define the nonprompt B model 

   RooDecay decay_nonprompt1("decay_nonprompt1", "decay", *time, tau0, effgaussmodel, RooDecay::SingleSided );       
   RooDecay decay_nonprompt2("decay_nonprompt2", "decay", *time, tau, effgaussmodel, RooDecay::SingleSided );
   RooRealVar fsl("fsl", "fraction of nonprompt component", 0.8, 0., 1.);
   RooAddPdf decay_nonprompt("nonp", "nonp", RooArgList( decay_nonprompt1,  decay_nonprompt2), fsl);


   RooRealVar tau1( "tau1", "tau1", 0.56, -1., 1.);
   RooRealVar tau2( "tau2", "tau2", 0.56, -1., 1.);

   RooDecay decay_wpv1("decay_wpv1", "decay", *time, tau1, tm, RooDecay::DoubleSided ); //convolv

   RooDecay decay_wpv2("decay_wpv2", "decay", *time, tau2, tm, RooDecay::DoubleSided ); //convolv

   RooRealVar fw("fwpv","fwpv",0.5, 0., 1.) ;
   RooAddPdf wpvPdf("wpvPdf", "wpvPdf",  RooArgList(decay_wpv1, decay_wpv2), fw );

    RooRealVar c1("c1", "c1", 4.66666e-02,  -10.0, 10.0); //org
    RooRealVar c2("c2", "c2", -3.54591e-01,  -10., 10.0); //org
//    RooAbsPdf * wpvPdf = new RooChebychev("wpvPdf","wpvPdf", *time, RooArgSet(c1,c2));

    RooDataSet data("data", "data", *time, Import(*tree));

    /// fractions ///

    RooRealVar f1("f1","f1",0.8, 0., 1.) ;
    RooRealVar f2("f2","f2",0.8, 0., 1.) ;

//    RooRealVar fprompt("fprompt","fprompt",0.8, 0., 1.) ;
//    RooRealVar fll("fll","fll",0.8, 0., 1.) ;

    RooFormulaVar fprompt("fprompt", "fprompt", "f1", RooArgList(f1));
    RooFormulaVar fll("fll", "fll", "(1-f2)*f1", RooArgList(f1, f2));
    RooFormulaVar fwpv("fwpv", "fwpv", "(1-f2)*(1-f1)", RooArgList(f1, f2));

    RooAddPdf model("model", "model", RooArgList(gaussmodel, decay_nonprompt), RooArgList(fprompt));
    model.fitTo(data, RooFit::Optimize(false));
    RooPlot *xframe = time->frame(100);
    data.plotOn(xframe);
    model.plotOn(xframe, Components(gaussmodel),  LineStyle(ELineStyle::kDashed), LineColor(kCyan));
    model.plotOn(xframe, Components(decay_nonprompt), LineStyle(ELineStyle::kDashed), LineColor(kRed));
 //   model.plotOn(xframe, Components(decay_tm),  LineStyle(ELineStyle::kDashed), LineColor(kGreen));
    model.plotOn(xframe, Components(wpvPdf), LineStyle(ELineStyle::kDashed), LineColor(kGreen));
//    model.paramOn(xframe, Format("NE",FixedPrecision(5)));
    model.plotOn(xframe);
  
 /*   decay_prompt.fitTo(data, RooFit::Optimize(false));
    RooPlot *xframe = time->frame(100);
    data.plotOn(xframe);
    decay_prompt.plotOn(xframe);*/

    TCanvas * c22 = new TCanvas();
    c22->SetLogy();
    xframe->Draw();
    xframe->SetTitle("");

    gPad->SetBottomMargin(0.36);

    RooHist *xpull = xframe->pullHist();
    RooPlot *Pullframe = time->frame();

    Pullframe->addPlotable(xpull, "P");

    Pullframe->GetXaxis()->SetTitle("treco");
    Pullframe->GetYaxis()->SetNdivisions(205, kTRUE);
    Pullframe->GetYaxis()->SetTickLength(0.083);
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 1.0 - 0.06, 0, 0, 0.);
    pad2->SetTopMargin(1.0 - 0.36);
    pad2->SetFillStyle(0);
    pad2->Draw();
    pad2->cd();
    Pullframe->Draw();
}

