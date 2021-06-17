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


void treco(std::string filename="h.root"){

  TFile *theFile = new TFile(filename.c_str());
  TTree *tree = (TTree *)theFile->Get("DecayTree");

  using namespace std;
  using namespace RooFit;  

  RooRealVar* time =new RooRealVar("time","t_{reco}[ps]", -0.2, 1.5) ;
  RooRealVar* trecoe=new RooRealVar("trecoe","terr[ps]",0.003,0.26);
  RooTruthModel* CM =new RooTruthModel("CM","",*time); // to be added!
  
  RooRealVar* meant=new RooRealVar("meant", "mean", 0, -0.01, 0.2);
  RooRealVar* sigmat=new RooRealVar("sigmat","sigma",0.05, -1., 1.);

  RooGaussian *gausspdf1 = new RooGaussian("gausspdf1", "gausspdf1",  *time, *meant , *sigmat ) ;
 // RooGaussian gausspdf("gausspdf1", "gausspdf1",  *time, *meant , *sigmat ) ;

  RooRealVar* sigmat2=new RooRealVar("sigmat2","sigma2",0.09, 0, 1.2);
  RooGaussian *gausspdf2 = new RooGaussian("gausspdf2", "gausspdf2",  *time, *meant , *sigmat2 ) ;


  RooRealVar* gmm=new RooRealVar("gmm", "mean", -0.02, -0.1, 0.2 );
  RooRealVar* gms=new RooRealVar("gms","sigma",0.01, -1., 1.);
  RooGaussModel gm("gm","gm", *time,*gmm, *gms);//This is the Gaussian resolution


  //Define the nonprompt B model
   RooRealVar a1("a1", "a1",  -0.5, -3., 3.); // well for single side
   RooExponential expoPdf("expoPdf", "expoPdf", *time, a1);


  //Define wPV events
   RooRealVar*  Tau =new RooRealVar("Tau","Tau", 0.1, -0.1, 10.);
 //   RooRealVar*  Tau =new RooRealVar("Tau","Tau", 1.2);
 //  RooDecay wpvPdf("roodecay","",*time, *Tau, *CM, RooDecay::DoubleSided); //Delta function convoluted
   RooDecay wpvPdf("roodecay","",*time,*Tau, gm, RooDecay::DoubleSided);
   
   RooDataSet data("data", "data", *time, Import(*tree));

   RooRealVar* ncore=new RooRealVar("ncore","prompt fraction", 20000, 0, tree->GetEntries());
   RooRealVar* Nnonp=new RooRealVar("Nnonp","nonpromt b ", 200, 0., tree->GetEntries());
   RooRealVar* nwPV=new RooRealVar("nwPV","NwPV", 10,  0, tree->GetEntries()); //0.005

   RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
   RooAddPdf gausspdf("sig","Signal",RooArgList(*gausspdf1, *gausspdf2), sig1frac) ;

//   RooAddPdf model("model", "model", RooArgList(gausspdf, expoPdf), bkgfrac);
  //  RooAddPdf model("model", "model", RooArgList(gausspdf, expoPdf), RooArgList(*ncore, *Nnonp));
// RooAddPdf model("model", "model", RooArgList(gausspdf, wpvPdf), RooArgList(*ncore, *nwPV));
    RooAddPdf model("model", "model", RooArgList(gausspdf, wpvPdf, expoPdf), RooArgList(*ncore, *nwPV, *Nnonp));
    model.fitTo(data, RooFit::Optimize(false));
    RooPlot *xframe = time->frame(100);
    data.plotOn(xframe);
    model.plotOn(xframe, Components(gausspdf), LineStyle(ELineStyle::kDashed), LineColor(kCyan));
    model.plotOn(xframe, Components(expoPdf),  LineStyle(ELineStyle::kDashed), LineColor(kGreen));
    model.plotOn(xframe, Components(wpvPdf), LineStyle(ELineStyle::kDashed), LineColor(kRed));
    model.paramOn(xframe, Format("NE",FixedPrecision(5)));
    model.plotOn(xframe);
  
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

