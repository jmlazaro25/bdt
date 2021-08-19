#include "EffPlotTools.hh"
#include "StyleTools.hh"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

#include <sstream>

using namespace EffPlotTools;
using namespace StyleTools;

double getMaximum(TH1D* h){
    double max = 0;
    for(int ibin = 1; ibin < h->GetNbinsX() + 1; ++ibin) {
        if(h->GetBinContent(ibin) > max) max = h->GetBinContent(ibin);
    }
    return max;
}

void drawRatio(
    TCanvas *c,
    vector<TH1D*> hists,
    vector<TString> labels,
    float legymax,
    float legymin,
    float rmax,
    float rmin,
    bool logy
) {

    c->cd();
    TPad *p1 = new TPad("p1", "p1", 0, 0.3, 1, 1);
    p1->SetLeftMargin(0.18);
    p1->SetTopMargin(0.10);
    p1->SetRightMargin(0.07);
    p1->SetBottomMargin(0.03);
    p1->Draw();
    if(logy) p1->SetLogy();
    p1->cd();
    TH1D* h1 = hists[0];
    double ymax = h1->GetMaximum();
    for(auto hist : hists) {
        if(hist->GetMaximum() > ymax) ymax = hist->GetMaximum();
            hist->SetLineWidth(3);
    }
    h1->GetXaxis()->SetLabelOffset(0.20);
    if(logy) h1->GetYaxis()->SetRangeUser(2, 1.5*ymax);
    else h1->GetYaxis()->SetRangeUser(0, 1.1*ymax);
    h1->Draw("axis");
    TLegend* leg = new TLegend(0.4, legymin, 0.6, legymax);
    for(unsigned int ih = 0; ih < hists.size(); ++ih) {
        if(logy) hists[ih]->GetYaxis()->SetRangeUser(2, 1.5*ymax);
        else hists[ih]->GetYaxis()->SetRangeUser(0, 1.1*ymax);
        hists[ih]->DrawCopy("histsame");
        leg->AddEntry(hists[ih], labels[ih], "L");
    }
    StyleTools::SetLegendStyle(leg);
    leg->Draw("same");
    c->cd();
    TPad *p2 = new TPad("p2", "p2", 0, 0, 1, 0.3);
    p2->SetLeftMargin(0.18);
    p2->SetTopMargin(0.00);
    p2->SetRightMargin(0.07);
    p2->SetBottomMargin(0.30);
    p2->SetGridy(1);
    p2->Draw();
    p2->cd();
    ymax = 0;
    for(unsigned int ih = 1; ih < hists.size(); ++ih) {
        TH1D* hist = hists[ih];
        hist->Divide(h1);
        if (getMaximum(hist) > ymax) ymax = getMaximum(hist);
        hist->SetTitleSize  (0.12, "Y");
        hist->SetTitleOffset(0.60, "Y");
        hist->SetTitleSize  (0.12, "X");
        hist->SetLabelSize  (0.10, "X");
        hist->SetLabelSize  (0.08, "Y");
        hist->GetYaxis()->SetTitleFont(62);
        hist->GetYaxis()->CenterTitle(kTRUE);
        hist->GetXaxis()->SetTitleFont(62);
        hist->GetYaxis()->SetNdivisions(305);
        hist->GetYaxis()->SetTitle("Ratio");
    }
    for(unsigned int ih = 1; ih < hists.size(); ++ih) {
        TH1D* hist = hists[ih];
        hist->GetYaxis()->SetRangeUser(rmin, rmax);
        hist->Draw("histsame");
    }
    double xmin = h1->GetXaxis()->GetXmin();
    double xmax = h1->GetXaxis()->GetXmax();
    TLine *l = new TLine(xmin, 1, xmax, 1);
    l->SetLineWidth(3);
    l->Draw("same");
    c->cd();
}

void SetupColors(){
    const unsigned num = 5;
    const int bands = 255;
    int colors[bands];
    double stops[num] = {0.00, 0.34, 0.61, 0.84, 1.00};
    double red[num] = {0.50, 0.50, 1.00, 1.00, 1.00};
    double green[num] = {0.50, 1.00, 1.00, 0.60, 0.50};
    double blue[num] = {1.00, 1.00, 0.50, 0.40, 0.50};
    int fi = TColor::CreateGradientColorTable(num, stops, red, green, blue, bands);
    for(int i = 0; i < bands; ++i){
        colors[i] = fi+i;
    }
    gStyle->SetNumberContours(bands);
    gStyle->SetPalette(bands, colors);
}

// Macro to plot variables stored in flat trees
void plot_bdt_recoilPt_0(
    const TString &bdtName1,
    const TString &bdtName2,
    const TString &flatDir1,
    const TString &flatDir2,
    const TString &outDir,
    const TString &withSelect,
    const double &bkgEff
) {

    // Background and signal names in the trees
    vector<TString> procs = {"0.001", "0.01", "0.1", "1.0", "bkg"};

    // Labels for the plots
    ostringstream ss;
    ss << bkgEff;
    TString effstr = ss.str();

    map<TString, TString> proclabel;
    proclabel["0.001"] = "m_{A'} = 0.001 GeV";
    proclabel["0.01"] = "m_{A'} = 0.01 GeV";
    proclabel["0.1"] = "m_{A'} = 0.1 GeV";
    proclabel["1.0"] = "m_{A'} = 1 GeV";
    proclabel["bkg"] = "bkg";
    proclabel["base"] = "All events";
    proclabel["fiducial"] = "Fiducial events";
    proclabel["nonfiducial"] = "Non-fiducial events";
    proclabel["trigger"] = "Trigger selected fiducial events";
    proclabel[bdtName1] = bdtName1 + " #varepsilon(bkg) = " + effstr;
    //proclabel[bdtName1] = "gabrielle (v9) #varepsilon(bkg) = " + effstr;
    proclabel[bdtName2] = bdtName2 + " #varepsilon(bkg) = " + effstr;
    //proclabel[bdtName2] = "gabrielle (v12) #varepsilon(bkg) = " + effstr;

    // Variable bin widths
    int n_bins = 14;
    int xbins[15] = {0, 2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100};
    Double_t xbins_dbl[n_bins + 1];
    for(size_t j = (n_bins + 1); j --;) {
        xbins_dbl[j] = xbins[j];
    }

    // Set up colors (defined in StyleTools)
    map<TString, unsigned int> colors;
    colors["bkg"] = color_comp1;
    colors["0.001"] = color_comp2;
    colors["0.01"] = color_comp3;
    colors["0.1"] = color_comp4;
    colors["1.0"] = color_comp5;
    colors["base"] = color_comp1;
    colors["fiducial"] = color_comp1;
    colors["nonfiducial"] = color_comp1;
    colors["trigger"] = color_comp1;
    colors[bdtName1] = color_comp2;
    colors[bdtName2] = color_comp3;

    // Define the selections we want to apply
    vector<TString> sels = {withSelect, bdtName1, bdtName2};

    // Cuts corresponding to different selections
    map<TString, TString> sel;
    sel["base"] = "1==1";
    sel["fiducial"] = sel["base"] + " && fiducial==1";
    sel["nonfiducial"] = sel["base"] + " && fiducial==0";
    sel["trigger"] = sel["base"] + " && trigPass==1";

    // Calculate cut values corresponding to the chosen selection
    for(auto seln : sels) {
        if(seln == bdtName1) {
            TH1D* bkghist = 0;

            // Read tree from input file
            TFile* file = new TFile(flatDir1 + "/bkg_" + bdtName1 + "_eval.root");
            //TFile* file = new TFile(flatDir1 + "/bkg_" + bdtName1 + "_eval_0.root");
            TTree* tree = (TTree*)file->Get("EcalVeto");
            assert(tree);

            // Set up histogram to be filled with disc values
            TH1D* dvhist = new TH1D(seln + "_bkg", "", 5000, 0.0, 1.0);

            // Read off disc values from tree into histogram
            tree->Draw("discValue_" + bdtName1 + ">>" + seln + "_bkg", sel[withSelect]);
            //tree->Draw("discValue_fernand>>" + seln + "_bkg", sel[withSelect]); //" + bdtName1 + ">>" + seln + "_bkg", sel[withSelect]);

            // Store background histogram
            bkghist = dvhist;

            // Get cut value for selected background efficiency
            float cutval = getCutValueForEfficiency(bkghist, bkgEff)[0];
            cout << "Cut value for BDT version (" + bdtName1 + "): cutval = " << cutval;
            sel[seln] = sel[withSelect] + " && discValue_" + bdtName1 + " > " + cutval;
            // sel[seln] = sel[withSelect] + " && discValue_fernand > " + cutval; //" + bdtName1 + " > " + cutval;
        }

        else if(seln == bdtName2) {
            TH1D* bkghist = 0;

            // Read tree from input file
            TFile* file = new TFile(flatDir2 + "/bkg_" + bdtName2 + "_eval.root");
            TTree* tree = (TTree*)file->Get("EcalVeto");
            assert(tree);

            // Set up histogram to be filled with disc values
            TH1D* dvhist = new TH1D(seln + "_bkg", "", 5000, 0.0, 1.0);

            // Read off disc values from tree into histogram
            tree->Draw("discValue_" + bdtName2 + ">>" + seln + "_bkg", sel[withSelect]);

            // Store background histogram
            bkghist = dvhist;

            // Get cut value for selected background efficiency
            float cutval = getCutValueForEfficiency(bkghist, bkgEff)[0];
            cout << "\nCut value for BDT version (" + bdtName2 + "): cutval = " << cutval;
            sel[seln] = sel[withSelect] + " && discValue_" + bdtName2 + " > " + cutval;
        }

        else
            continue;
    }

    // Make nice looking plots
    SetTDRStyle();
    SetupColors();

    gStyle->SetLabelSize(0.03, "XYZ");

    // Create the output directory
    gSystem->mkdir(outDir, true);

    // Loop over processes
    for(auto proc : procs) {
        cout << "\nGetting data for process with label (" + proc + ")..." << endl;

        // Vector to hold all the histograms
        vector<TH1D*> hists;
        double ymax = -1.;
      
        // Loop over selections
        for(auto seln : sels) {

            // Read trees from input file
            TFile* file = NULL;
            TTree* tree = NULL;

            if(seln == bdtName1) {
                cout << "Getting data for BDT version (" + bdtName1 + ")..." << endl;
                file = new TFile(flatDir1 + "/" + proc + "_" + bdtName1 + "_eval.root");
                //file = new TFile(flatDir1 + "/" + proc + "_" + bdtName1 + "_eval_0.root");
                tree = (TTree*)file->Get("EcalVeto");
            }

            else {
                cout << "Getting data for BDT version (" + bdtName2 + ")..." << endl;
                file = new TFile(flatDir2 + "/" + proc + "_" + bdtName2 + "_eval.root");
                tree = (TTree*)file->Get("EcalVeto");
            }

            assert(tree);

            // Set up histogram with the correct binning
            TH1D* hist = new TH1D("h_recoilPt_" + proc + "_" + seln, "", 100, 0, 100);
	    hist->SetBins(n_bins, xbins_dbl);
		
            hist->SetLineColor(colors[seln]);

            // Draw the plot using the given data
            tree->Draw("recoilPT >> h_recoilPt_" + proc + "_" + seln, sel[seln], "hist");

            // Add overflow contents to the last bin
            addOverFlow(hist);

            // Set up axis labels and colors
            hist->GetXaxis()->SetTitle("p_{T}(recoil e^{-}) (MeV)");
            hist->GetXaxis()->SetTitleSize(0.04);
            hist->GetYaxis()->SetTitle("Entries");
            hist->GetYaxis()->SetTitleSize(0.05);
            hist->SetLineWidth(3);
            hist->SetLineColor(colors[seln]);

            // Store maximum y value
            if(hist->GetMaximum() > ymax) ymax = hist->GetMaximum();

            // Add it to the list
            hists.push_back(hist);
        }

        // Set up the canvas and legend
        TCanvas* c = MakeCanvas("c", "", 600, 600);
        TLegend* leg = new TLegend(0.45, 0.6, 0.85, 0.85);
        leg = new TLegend(0.2, 0.6, 0.6, 0.85);
        SetLegendStyle(leg);
        c->SetLogy();

        vector<TString> labels;

        // Draw histograms for all selections on the canvas
        for(auto* hist : hists) {

            // Set y-axis scales
            hist->GetYaxis()->SetRangeUser(1, 1.5*ymax);

            // Draw it!
            hist->Draw("histsame");

            // Add this process to the legend
            TString hname(hist->GetName());
            TString procname = TString(hname(hname.Last('_') + 1, hname.Length()));
            labels.push_back(proclabel[procname]);
        }

        // Draw the ratio plots
        if(proc == "0.001")
            drawRatio(c, hists, labels, 0.85, 0.65, 1.0, 0.3, true);
        else if(proc == "0.01")
            drawRatio(c, hists, labels, 0.85, 0.65, 1.0, 0.3, true);
        else if(proc == "0.1")
            drawRatio(c, hists, labels, 0.3, 0.1, 1.0, 0.3, true);
        else if(proc == "1.0")
            drawRatio(c, hists, labels, 0.3, 0.1, 1.0, 0.9, true);
        else
            drawRatio(c, hists, labels, 0.85, 0.65, 0.002, 0.00005, true);

        // Draw the legend
        leg->Draw("same");

        // Draw LDMX text
        LDMX_lumi(c, 0);

        // Save the plot
        c->SaveAs(outDir + "/" + withSelect + "_" + bdtName1 + "_" + bdtName2 + "_" + proc + "_recoilPt.pdf");
    }
}
