//plot a single ROC curve set (init bdt analysis)

#include "StyleTools.hh"
#include "EffPlotTools.hh"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TMarker.h"

using namespace StyleTools;
using namespace EffPlotTools;


// Macro to plot ROC curves
// To run: root -l -q -b plotROC.C+ 
void plot1ROC(const TString inputdir="/home/jmlazaro/research/bdt/test_bdt/evals", // input directory with flat trees
             const TString outputdir="/home/jmlazaro/research/bdt/test_bdt/plots", // output directory to save plots into
             const TString groupings="m",  
             const TString tree_name="gabrielle_back_v1",  
             const bool zoom = true, // create a zoomed plot (ranges defined below)
             const bool log = false // create a zoomed plot (ranges defined below)
             )
{

  const float intcut = 0.90;
  const float inteff = 0.001;

  Double_t x1, y1, x2, y2;
  x1 = 0.5;
  y1 = 0.2;
  x2 = 0.9;
  y2 = 0.5;

  TString bkg, sig1, sig2, sig3, sig4;
  map<TString,TString> proclabel;

  if(groupings=="m") {
      // Bkg + signal names in the trees
      bkg = "pn";
      sig1 = "0.001";
      sig2 = "0.01";
      sig3 = "0.1";
      sig4 = "1.0";

      // Hadd the input flat trees if they haven't been already
      proclabel[bkg] = "Photo-nuclear";
      proclabel[sig1] = "m_{A'} = 0.001 GeV";
      proclabel[sig2] = "m_{A'} = 0.01 GeV";
      proclabel[sig3] = "m_{A'} = 0.1 GeV";
      proclabel[sig4] = "m_{A'} = 1.0 GeV";
    }

  if(groupings=="e") {
      // Bkg + signal names in the trees
      bkg = "neutron";
      sig1 = "e0.2";
      sig2 = "e0.5";
      sig3 = "e1";
      sig4 = "e2";

      // Hadd the input flat trees if they haven't been already
      proclabel[bkg] = "Neutrons";
      proclabel[sig1] = "E_{#gamma} = 0.2 GeV";
      proclabel[sig2] = "E_{#gamma} = 0.5 GeV";
      proclabel[sig3] = "E_{#gamma} = 1.0 GeV";
      proclabel[sig4] = "E_{#gamma} = 2.0 GeV";
    }

  if(groupings=="t") {
      // Bkg + signal names in the trees
      bkg = "neutron";
      sig1 = "t0";
      sig2 = "t15";
      sig3 = "t25";
      sig4 = "t55";

      // Hadd the input flat trees if they haven't been already
      proclabel[bkg] = "Neutrons";
      proclabel[sig1] = "#theta_{#gamma} = 0#circ";
      proclabel[sig2] = "#theta_{#gamma} = 15#circ";
      proclabel[sig3] = "#theta_{#gamma} = 25#circ";
      proclabel[sig4] = "#theta_{#gamma} = 55#circ";
    }

  vector<TString> procs = {bkg,sig1,sig2,sig3,sig4};

  // Set up colors (defined in StyleTools)
  map<TString,unsigned int> colors;
  colors[bkg] = color_comp1;
  colors[sig1] = color_comp2;
  colors[sig2] = color_comp3;
  colors[sig3] = color_comp4;
  colors[sig4] = color_comp5;

  // Cuts corresponding to different selections
  map<TString,TString> sel;
  sel["base"] = "1==1";
  sel["fiducial"] = sel["base"] + " && fiducial==1";
  sel["nonfiducial"] = sel["base"] + " && fiducial==0";
  sel["trig"] = sel["base"];
  sel["sr"] = sel["trig"] + " && sqrt(recoilPx*recoilPx + recoilPy*recoilPy + recoilPz*recoilPz) < 1200";

  //sel["trig"] = sel["fiducial"] + " && trigPass == 1";
  //sel["base"] = "pnWeight>0.2";
  //sel["sr"] = sel["fiducial"] + " && sqrt(recoilPx*recoilPx + recoilPy*recoilPy + recoilPz*recoilPz) < 1200";


  // Define the selections we want to apply
  vector<TString> sels = {"base"};


  // Make nice looking plots
  SetTDRStyle();

  gSystem->mkdir(outputdir, true);

  // Loop over selections
  for(auto seln : sels) {
    TH1D* bkghist = 0;
    vector<TH1D*> sighist;

    // Loop over processes
    for(auto proc : procs) {
      // Read tree from input files
      TFile* file = new TFile(inputdir+"/"+proc+"_eval.root");
      TTree* tree = (TTree*)file->Get(tree_name+"_Veto");
      assert(tree);

      // Set up histogram to be filled with the disc values
      TH1D* hist = new TH1D("discvalue_"+seln+"_"+proc,"",5000,0.0,1.0);

      // Read off disc values from tree into the histogram
      tree->Draw("discValue_"+tree_name+">>discvalue_"+seln+"_"+proc, sel[seln]);

      // Store histogram as background or signal histogram
      if(proc==bkg) {
        bkghist = hist;
      }
      else {
        sighist.push_back(hist); //vector of histograms so need to push back
      }
    }


    // Get background eff for intcut
    cout << "neutron" << endl;
    cout << "intcut = " << intcut << endl;
    float bkgeff = getEfficiencyForCutValue(bkghist, intcut)[0];
    cout << "intcut: eff(bkg) = " << bkgeff << endl;

    // Get cut value for inteff
    cout << "inteff = " << inteff << endl;
    float cutval = getCutValueForEfficiency(bkghist, inteff)[0];
    cout << "cut for inteff = " << cutval << endl;
    cout << " " << endl;

    // Set up the canvas and legend
    TCanvas* c1 = MakeCanvas("c1","",600,600);
    TLegend* leg = new TLegend(x1,y1,x2,y2);
    //TLegend* leg = new TLegend(0.45,0.2,0.8,0.65);
    leg->SetTextSize(0.02);
    SetLegendStyle(leg);


    // Loop over signals and compute ROC curve for each
    int isig = 0;
    for(auto* sighist : sighist) {
      TString hname(sighist->GetName());

      /*
      if (isig == 0) {
            TString procname = TString('u');
        } else if (time < 20) {
            TString procname = TString('_');
        } else {
            TString procname = TString(hname(hname.Last('_')+1,hname.Length()));
        }
      */
      TString procname = TString(procs[isig + 1]);

      // Make the ROC curve
      TGraph* rocgr = computeROCCurve(sighist,bkghist,"",false,false,true);

      // Get signal eff for intcut
      float sigeff = getEfficiencyForCutValue(sighist, intcut)[0];
      cout << "eff (" << procname << ") for intcut = " << sigeff << endl;
      cout << "eff (" << procname << ") for " << cutval << " cut = " << getEfficiencyForCutValue(sighist, cutval)[0] << endl;
      cout << " " << endl;

      // Set up colors and style
      rocgr->SetLineColor(colors[procname]);
      rocgr->SetMarkerColor(colors[procname]);
      rocgr->SetLineWidth(3);
      //rocgr->SetLineStyle(2); //makes dashed line
      
      // Add this signal to the legend
      leg->AddEntry(rocgr,proclabel[procname],"L");  // + ", HCal BDT","L");
      //leg->AddEntry(rocgr,proclabel[procname] + ", bdt_2","L");
      //leg->AddEntry(rocgr,proclabel[procname],"L");  // + ", HCal BDT","L");

      // Set axis titles
      rocgr->GetXaxis()->SetTitle("#varepsilon(background)");
      rocgr->GetYaxis()->SetTitle("#varepsilon(signal)");


      // Zoom in on and set log for the x axis (bkg eff) if desired
      if(zoom) rocgr->GetXaxis()->SetLimits(0.0,0.1);
      if(zoom) rocgr->GetYaxis()->SetRangeUser(0.3,1);
      if(zoom) rocgr->GetXaxis()->SetLimits(0.0,0.02);
      if(zoom) rocgr->GetYaxis()->SetRangeUser(0.4,1);
      if(zoom) rocgr->GetXaxis()->SetLimits(0.0,0.1);
      if(zoom) rocgr->GetYaxis()->SetRangeUser(0.3,1);
      if(log) {
        rocgr->GetXaxis()->SetLimits(0.0001,0.1);
        c1->SetLogx();
      }


      // If it's the first process, draw the axis as well
      c1->cd();
      if(isig==0) rocgr->Draw("AC"); //A to draw axis, C to draw smooth curve
      else rocgr->Draw("Csame"); //same to put on same canvas


      // Marker showing the position of the 0.99 cut on the curve
      TMarker* m = new TMarker(bkgeff,sigeff,29); //22 for filled in triangle
      m->SetMarkerColor(6);        //(colors[procname]);
      m->SetMarkerSize(1.7);
      m->Draw("same");

      isig++; //important to iterate through different mass points
    }
    cout << endl;

    // Draw the legend
    c1->cd();
    leg->Draw("same");
    // Draw LDMX text
    LDMX_lumi(c1,0);
    // Save the plot
    c1->SaveAs(outputdir+"/" + seln + (zoom ? "_zoom" : "") + (log ? "_log" : "") + "_ROC.pdf");
    c1->Close();
  }
}
