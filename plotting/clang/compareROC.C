//plot two ROC curve sets (comparing bdts)

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
void compareROC(const TString inputdirnom="/nfs/slac/g/ldmx/users/jmlazaro/bdt/gabrielle_bdt/evals",
	            const TString inputdirnew="/nfs/slac/g/ldmx/users/jmlazaro/bdt/back_v1_bdt/evals",
	            const TString outputdirnom="/nfs/slac/g/ldmx/users/jmlazaro/bdt/gabrielle_bdt/plots",
	            const TString outputdirnew="/nfs/slac/g/ldmx/users/jmlazaro/bdt/back_v1_bdt/plots",
                const TString tree_name_nom="gabrielle",
                const TString tree_name_new="back_v1", 
                const bool zoom = true, // create a zoomed plot (ranges defined below)
                const bool log = false // create a zoomed plot (ranges defined below)
                )
{

    const TString outputdir = inputdirnew;  // output directory to save plots into

  const TString nomStr = "Gabrielle";
  const TString newStr = "Backv1";

  const float intcut = 0.90;
  const float inteff = 0.001;

  // Bkg + signal names in the trees
  vector<TString> procs = {"pn", "0.001", "0.01", "0.1", "1.0"};

  // Hadd the input flat trees if they haven't been already
  map<TString,TString> proclabel;
  proclabel["pn"    ] = "Photo-nuclear";
  proclabel["0.001" ] = "m_{A'} = 0.001 GeV";
  proclabel["0.01"  ] = "m_{A'} = 0.01 GeV";
  proclabel["0.1"   ] = "m_{A'} = 0.1 GeV";
  proclabel["1.0"   ] = "m_{A'} = 1.0 GeV";

  // Set up colors (defined in StyleTools)
  map<TString,unsigned int> colors;
  colors["pn"    ] = color_comp1;
  colors["0.001" ] = color_comp2;
  colors["0.01"  ] = color_comp3;
  colors["0.1"   ] = color_comp4;
  colors["1.0"   ] = color_comp5;
  /*
  // Hadd the input flat trees if they haven't been already
  map<TString,TString> proclabel;
  proclabel["neutron"] = "Neutrons";
  proclabel["e0.2"] = "E_{#gamma} = 0.2 GeV";
  proclabel["e0.5"] = "E_{#gamma} = 0.5 GeV";
  proclabel["e1"] = "E_{#gamma} = 1.0 GeV";
  proclabel["e2"] = "E_{#gamma} = 2.0 GeV";

  // Set up colors (defined in StyleTools)
  map<TString,unsigned int> colors;
  colors["neutron"] = color_comp1;
  colors["e0.2"] = color_comp2;
  colors["e0.5"] = color_comp3;
  colors["e1"] = color_comp4;
  colors["e2"] = color_comp5;
  */

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

  gSystem->mkdir(outputdirnom, true);
  gSystem->mkdir(outputdirnew, true);


  // Loop over selections
  for(auto seln : sels) {
    TH1D* bkghistnom = 0;
    TH1D* bkghistnew = 0;
    vector<TH1D*> sighistsnom;
    vector<TH1D*> sighistsnew;

    // Loop over processes
    for(auto proc : procs) {
      // Read tree from input files
      TFile* filenom = new TFile(inputdirnom+"/"+proc+"_eval.root");
      TTree* treenom = (TTree*)filenom->Get(tree_name_nom+"_Veto");
      assert(treenom);

      TFile* filenew = new TFile(inputdirnew+"/"+proc+"_eval.root");
      TTree* treenew = (TTree*)filenew->Get(tree_name_new+"_Veto");
      assert(treenew);


      // Set up histogram to be filled with the disc values
      TH1D* histnom = new TH1D("discvalue_nom_"+seln+"_"+proc,"",5000,0.0,1.0);
      TH1D* histnew = new TH1D("discvalue_new_"+seln+"_"+proc,"",5000,0.0,1.0);


      // Read off disc values from tree into the histogram
      treenom->Draw("discValue_"+tree_name_nom+">>discvalue_nom_"+seln+"_"+proc,sel[seln]);
      treenew->Draw("discValue_"+tree_name_new+">>discvalue_new_"+seln+"_"+proc,sel[seln]);


      // Store histogram as background or signal histogram
      if(proc==procs[0]) {
        bkghistnom = histnom;
        bkghistnew = histnew;
      }
      else {
        sighistsnom.push_back(histnom); //vector of histograms so need to push back
        sighistsnew.push_back(histnew);
      }
    }


    // Get background eff for intcut
    cout << "neutron" << endl;
    cout << "intcut = " << intcut << endl;
    float bkgeffnom = getEfficiencyForCutValue(bkghistnom, intcut)[0];
    cout << "intcut: eff(bkg) nom = " << bkgeffnom << endl;
    float bkgeffnew = getEfficiencyForCutValue(bkghistnew, intcut)[0];
    cout << "intcut: eff(bkg) new = " << bkgeffnew << endl;


    // Get cut value for 1e-4 bkg eff
    cout << "inteff = " << inteff << endl;
    float cutval = getCutValueForEfficiency(bkghistnom, inteff)[0];
    cout << "cut for inteff = " << cutval << endl;
    cout << " " << endl;

    // Set up the canvas and legend
    TCanvas* c1 = MakeCanvas("c1","",600,600);
    TLegend* leg = new TLegend(0.4,0.15,0.8,0.6);
    leg->SetTextSize(1);
    SetLegendStyle(leg);


    // Loop over signals and compute ROC curve for each
    int isig = 0;
    for(auto* sighistnom : sighistsnom) {
      TH1D* sighistnew = sighistsnew[isig];
      TString hname(sighistnom->GetName());
      TString procname = TString(hname(hname.Last('_')+1,hname.Length()));

      // Make the ROC curve
      TGraph* rocgrnom = computeROCCurve(sighistnom,bkghistnom,"",false,false,true);
      TGraph* rocgrnew = computeROCCurve(sighistnew,bkghistnew,"",false,false,true);

      // Get signal eff for intcut
      float sigeffnom = getEfficiencyForCutValue(sighistnom, intcut)[0];
      float sigeffnew = getEfficiencyForCutValue(sighistnew, intcut)[0];
      cout << "eff nom (" << procname << ") for intcut = " << sigeffnom << endl; 
      cout << "eff new (" << procname << ") for intcut = " << sigeffnew << endl; 

      cout << "eff nom (" << procname << ") for " << cutval << " cut = " << getEfficiencyForCutValue(sighistnom, cutval)[0] << endl;
      cout << "eff new (" << procname << ") for " << cutval << " cut = " << getEfficiencyForCutValue(sighistnew, cutval)[0] << endl;
      cout << " " << endl;

      // Set up colors and style
      rocgrnom->SetLineColor(colors[procname]);
      rocgrnom->SetMarkerColor(colors[procname]);
      rocgrnom->SetLineWidth(3);
      rocgrnom->SetLineStyle(2); //makes dashed line
      
      rocgrnew->SetLineColor(colors[procname]);
      rocgrnew->SetMarkerColor(colors[procname]);
      rocgrnew->SetLineWidth(3);
      //rocgrnew->SetLineStyle(2);
 

      // Add this signal to the legend
      leg->AddEntry(rocgrnom,proclabel[procname] + ", "+nomStr,"L");
      leg->AddEntry(rocgrnew,proclabel[procname] + ", "+newStr,"L");


      // Set axis titles
      rocgrnom->GetXaxis()->SetTitle("#varepsilon(background)");
      rocgrnom->GetYaxis()->SetTitle("#varepsilon(signal)");


      // Zoom in on and set log for the x axis (bkg eff) if desired
      if(zoom) rocgrnom->GetXaxis()->SetLimits(0.0,0.1);
      if(zoom) rocgrnom->GetYaxis()->SetRangeUser(0.3,1);
      if(zoom) rocgrnew->GetXaxis()->SetLimits(0.0,0.1);
      if(zoom) rocgrnew->GetYaxis()->SetRangeUser(0.3,1);
      if(log) {
        rocgrnom->GetXaxis()->SetLimits(0.0001,0.1);
        c1->SetLogx();
      }


      // If it's the first process, draw the axis as well
      c1->cd();
      if(isig==0)
	rocgrnom->Draw("AC"); //A to draw axis, C to draw smooth curve
      else rocgrnom->Draw("Csame"); //same to put on same canvas
      rocgrnew->Draw("Csame");


      // Marker showing the position of the 0.99 cut on the curve
      TMarker* mnom = new TMarker(bkgeffnom,sigeffnom,22); //22 for filled in triangle
      mnom->SetMarkerColor(colors[procname]);
      mnom->SetMarkerSize(1.2);
      mnom->Draw("same");

      TMarker* mnew = new TMarker(bkgeffnew,sigeffnew,29); //29 for filled in star
      mnew->SetMarkerColor(colors[procname]);
      mnew->SetMarkerSize(1.3);
      mnew->Draw("same");

      isig++; //important to iterate through different mass points
    }
    cout << endl;

    // Draw the legend
    c1->cd();
    leg->Draw("same");
    // Draw LDMX text
    LDMX_lumi(c1,0);
    // Save the plot
    c1->SaveAs(outputdirnom+"/" + seln + (zoom ? "_zoom" : "") + (log ? "_log" : "") + "_compareROC_"+tree_name_nom+"v"+tree_name_new+".pdf");
    c1->SaveAs(outputdirnew+"/" + seln + (zoom ? "_zoom" : "") + (log ? "_log" : "") + "_compareROC_"+tree_name_nom+"v"+tree_name_new+".pdf");
    c1->Close();
  }
}
