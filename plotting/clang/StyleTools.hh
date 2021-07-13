//--------------------------------------------------------------------------------------------------
//
// Class to implement a few simple style functions to make slightly better looking ROOT plots.
//
//--------------------------------------------------------------------------------------------------
#ifndef STYLETOOLS_HH
#define STYLETOOLS_HH

#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TASImage.h>
#include <iostream>
#include "map"


namespace StyleTools {

  enum Colors {color_data=1,                 // kBlack
               color_qcd=811,                // kSpring-9
               color_ttbar=866,              // kAzure+6
               color_ttbar1l=866,            // kAzure+6
               color_ttbar2l=625,            // kRed-7
               color_ttW=606,                // kMagenta-10
               color_ttZ=797,                // kOrange-3
               color_tW=876,                 // kViolet-4
               color_wjets=391,              // kYellow-9
               color_znunu=623,              // kRed-9
               color_other=391,              // kYellow-9
               color_T2tt=634,               // kRed+2
               color_T2tt_850_100=634,       // kRed+2
               color_T2tt_650_325=633,       // kRed+1
               color_T2tt_500_325=635,       // kRed+3
               color_T2tt_425_325=622,       // kRed-10
               color_T1tttt=602,             // kBlue+2
               color_T1tttt_1500_100=602,    // kBlue+2
               color_T1tttt_1200_800=601,    // kBlue+1
               color_comp1=634,              // kRed+2
               color_comp2=862,              // kAzure+2
               color_comp3=797,              // kOrange-3
               color_comp4=882,              // kViolet+2
               color_comp5=417,              // kGreen+3-3
               color_comp6=603,              // kBlue+3
               color_comp7=802,              // kOrange+2
               color_comp8=616,              // kMagenta
               color_comp9=600,              // kBlue
               color_comp10=434,             // kCyan+2
               color_comp11=800,             // kOrange
               color_comp12=417,             // kGreen+1
               color_comp13=632              // kRed
  };

  typedef std::map<TString, unsigned int> ColorMap;

  TCanvas*  MakeCanvas     (const char* name, const char *title, int dX = 500, int dY = 500);
  void      InitSubPad     (TPad* pad, int i);
  void      InitHist       (TH1 *hist, const char *xtit, const char *ytit  = "Number of Entries", int color = kBlack, int style = 0);
  void      InitGraph      (TGraph *gr, const char *title, const char *xtit, const char *ytit, double xmin, double xmax, double ymin, double ymax, int color = kBlack);
  void      SetLegendStyle (TLegend* leg);
  void      SetStyle       ();
  void      SetTDRStyle    ();
  void      CMS_lumi       (TPad* pad, int iPeriod=5, int iPosX=10, TString extraText = "Simulation", TString lumi_pass = "");
  void      LDMX_lumi      (TPad* pad, int iPosX=10, TString extraText = "Simulation");
  ColorMap  DefaultColors  ();
  
  TCanvas* newAlbum(Int_t numItems, const Char_t* name, const Char_t* title, Int_t width, Int_t height, Int_t rows, Int_t columns);
  TCanvas* drawAll(TObjArray* plots, TString name);
  int colorGetter(int counter); //For when you need quick new colors!
  TH1D* addOverFlow(TH1F* h);


TCanvas* MakeCanvas(const char* name, const char *title, int dX, int dY)
{
  // Start with a canvas
  TCanvas *canvas = new TCanvas(name,title,0,0,dX,dY);

  // General overall stuff
  canvas->SetFillColor      (0);
  canvas->SetBorderMode     (0);
  canvas->SetBorderSize     (10);

  // Set margins to reasonable defaults
  canvas->SetLeftMargin     (0.18);
  canvas->SetRightMargin    (0.05);
  canvas->SetTopMargin      (0.08);
  canvas->SetBottomMargin   (0.15);
  
  // Setup a frame which makes sense
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameLineWidth (3);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);

  return canvas;
}

void InitSubPad(TPad* pad, int i)
{
  pad->cd(i);
  TPad *tmpPad = (TPad*) pad->GetPad(i);
  tmpPad->SetLeftMargin  (0.18);
  tmpPad->SetTopMargin   (0.05);
  tmpPad->SetRightMargin (0.07);
  tmpPad->SetBottomMargin(0.15);
  return;
}

void InitHist(TH1 *hist, const char *xtit, const char *ytit, int color, int style)
{
  hist->SetXTitle(xtit);
  hist->SetYTitle(ytit);
  hist->SetLineColor(kBlack);
  hist->SetLineWidth(    3.);
  hist->SetFillColor(color );
  hist->SetFillStyle(style );
  hist->SetTitleSize  (0.055,"Y");
  hist->SetTitleOffset(1.600,"Y");
  hist->SetLabelOffset(0.014,"Y");
  hist->SetLabelSize  (0.040,"Y");
  hist->SetLabelFont  (62   ,"Y");
  hist->SetTitleSize  (0.055,"X");
  hist->SetTitleOffset(1.300,"X");
  hist->SetLabelOffset(0.014,"X");
  hist->SetLabelSize  (0.040,"X");
  hist->SetLabelFont  (62   ,"X");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize (1.3);
  hist->GetYaxis()->SetTitleFont(62);
  hist->GetXaxis()->SetTitleFont(62);
  hist->SetTitle("");  
  return;
}

void InitGraph(TGraph *gr, const char *title, const char *xtit, const char *ytit, double xmin, double xmax, double ymin, double ymax, int color)
{
  gr->SetTitle(title);
  gr->GetXaxis()->SetTitle(xtit);
  gr->GetYaxis()->SetTitle(ytit);
  gr->SetLineColor(color);
  gr->SetLineWidth(    3.);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(color);
  gr->SetMarkerSize (1.3);
  gr->GetYaxis()->SetTitleSize  (0.055);
  gr->GetYaxis()->SetTitleOffset(1.600);
  gr->GetYaxis()->SetLabelOffset(0.014);
  gr->GetYaxis()->SetLabelSize  (0.040);
  gr->GetYaxis()->SetLabelFont  (62   );
  gr->GetXaxis()->SetTitleSize  (0.055);
  gr->GetXaxis()->SetTitleOffset(1.300);
  gr->GetXaxis()->SetLabelOffset(0.014);
  gr->GetXaxis()->SetLabelSize  (0.040);
  gr->GetXaxis()->SetLabelFont  (62   );
  gr->GetYaxis()->SetTitleFont(62);
  gr->GetXaxis()->SetTitleFont(62);
  gr->GetXaxis()->SetLimits(xmin, xmax);
  gr->GetYaxis()->SetLimits(ymin, ymax);
  return;
}

void SetLegendStyle(TLegend* leg)
{
  leg->SetFillStyle (0);
  leg->SetFillColor (0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
}

void SetStyle()
{
  TStyle *MyStyle = new TStyle("New Style","Better than the default style :-)");
  gStyle = MyStyle;

  // Canvas
  MyStyle->SetCanvasColor     (0);
  MyStyle->SetCanvasBorderSize(10);
  MyStyle->SetCanvasBorderMode(0);
  MyStyle->SetCanvasDefH      (700);
  MyStyle->SetCanvasDefW      (700);
  MyStyle->SetCanvasDefX      (100);
  MyStyle->SetCanvasDefY      (100);

  // color palette for 2D temperature plots
  MyStyle->SetPalette(1,0);

  // Pads
  MyStyle->SetPadColor       (0);
  MyStyle->SetPadBorderSize  (10);
  MyStyle->SetPadBorderMode  (0);
  MyStyle->SetPadBottomMargin(0.13);
  MyStyle->SetPadTopMargin   (0.08);
  MyStyle->SetPadLeftMargin  (0.15);
  MyStyle->SetPadRightMargin (0.05);
  MyStyle->SetPadGridX       (0);
  MyStyle->SetPadGridY       (0);
  MyStyle->SetPadTickX       (1);
  MyStyle->SetPadTickY       (1);

  // Frames
  MyStyle->SetLineWidth(3);
  MyStyle->SetFrameFillStyle ( 0);
  MyStyle->SetFrameFillColor ( 0);
  MyStyle->SetFrameLineColor ( 1);
  MyStyle->SetFrameLineStyle ( 0);
  MyStyle->SetFrameLineWidth ( 3);
  MyStyle->SetFrameBorderSize(10);
  MyStyle->SetFrameBorderMode( 0);

  // Histograms
  MyStyle->SetHistFillColor(2);
  MyStyle->SetHistFillStyle(0);
  MyStyle->SetHistLineColor(1);
  MyStyle->SetHistLineStyle(0);
  MyStyle->SetHistLineWidth(3);
  MyStyle->SetNdivisions(505, "X");

  // Functions
  MyStyle->SetFuncColor(1);
  MyStyle->SetFuncStyle(0);
  MyStyle->SetFuncWidth(2);

  // Various
  MyStyle->SetMarkerStyle(20);
  MyStyle->SetMarkerColor(kBlack);
  MyStyle->SetMarkerSize (1.4);

  MyStyle->SetTitleBorderSize(0);
  MyStyle->SetTitleFillColor (0);
  MyStyle->SetTitleX         (0.2);

  MyStyle->SetTitleSize  (0.055,"X");
  MyStyle->SetTitleOffset(1.200,"X");
  MyStyle->SetLabelOffset(0.005,"X");
  MyStyle->SetLabelSize  (0.040,"X");
  MyStyle->SetLabelFont  (62   ,"X");

  MyStyle->SetStripDecimals(kFALSE);
  MyStyle->SetLineStyleString(11,"20 10");

  MyStyle->SetTitleSize  (0.055,"Y");
  MyStyle->SetTitleOffset(1.600,"Y");
  MyStyle->SetLabelOffset(0.010,"Y");
  MyStyle->SetLabelSize  (0.040,"Y");
  MyStyle->SetLabelFont  (62   ,"Y");

  MyStyle->SetTextSize   (0.055);
  MyStyle->SetTextFont   (62);

  MyStyle->SetStatFont   (62);
  MyStyle->SetTitleFont  (62);
  MyStyle->SetTitleFont  (62,"X");
  MyStyle->SetTitleFont  (62,"Y");

  MyStyle->SetOptStat    (0);

  return;
}

void SetTDRStyle()
{

  // Copied from https://ghm.web.cern.ch/ghm/plots/MacroExample/tdrstyle.C

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  // tdrStyle->SetErrorX(0.);

  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

// For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:
  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:
  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:
  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  //tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetNdivisions(505, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Legend font
  tdrStyle->SetLegendFont(42);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);
  //
  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->SetLineStyleString(11,"20 10");

  tdrStyle->cd();

  return;
}

void CMS_lumi(TPad* pad, int iPeriod, int iPosX, TString extraText, TString lumi_pass)
{

  // Global variables

  TString cmsText     = "CMS";
  float cmsTextFont   = 61;  // default is helvetic-bold

  bool writeExtraText = true;
  float extraTextFont = 52;  // default is helvetica-italics
  
  // text sizes and text offsets with respect to the top frame
  // in unit of the top margin size
  float lumiTextSize     = 0.6;
  float lumiTextOffset   = 0.2;
  float cmsTextSize      = 0.75;
  float cmsTextOffset    = 0.11; //0.1;  // only used in outOfFrame version

  float relPosX    = 0.045;
  float relPosY    = 0.035;
  float relExtraDY = 1.2;

  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;

  TString lumi_13TeV   = "2.3 fb^{-1}";
  TString lumi_13TeV_2 = "12.9 fb^{-1}";
  TString lumi_8TeV    = "19.7 fb^{-1}";
  TString lumi_7TeV    = "5.1 fb^{-1}";
  TString lumi_sqrtS   = "13 TeV";

  bool drawLogo      = false;

  bool outOfFrame    = false;
  if( iPosX/10==0 ) 
    {
      outOfFrame = true;
    }
  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  if( iPosX == 0  ) relPosX = cmsTextOffset;
  int align_ = 10*alignX_ + alignY_;

  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  //  float e = 0.025;

  pad->cd();

  TString lumiText;
  if( iPeriod==-1 ) {
      lumiText += lumi_pass;
      lumiText += " (13 TeV)";
  }
  else if( iPeriod==1 ) {
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
  }
  else if ( iPeriod==2 ) {
    lumiText += lumi_8TeV;
    lumiText += " (8 TeV)";
  }
  else if( iPeriod==3 ) {
    lumiText = lumi_8TeV; 
    lumiText += " (8 TeV)";
    lumiText += " + ";
    lumiText += lumi_7TeV;
    lumiText += " (7 TeV)";
  }
  else if ( iPeriod==4 ) {
    lumiText += lumi_13TeV;
    lumiText += " (13 TeV)";
  }
  else if ( iPeriod==5 ) {
    lumiText += lumi_13TeV_2;
    lumiText += " (13 TeV)";
  }
  else if ( iPeriod==7 ) { 
    if( outOfFrame ) lumiText += "#scale[0.85]{";
    lumiText += lumi_13TeV; 
    lumiText += " (13 TeV)";
    lumiText += " + ";
    lumiText += lumi_8TeV; 
    lumiText += " (8 TeV)";
    lumiText += " + ";
    lumiText += lumi_7TeV;
    lumiText += " (7 TeV)";
    if( outOfFrame) lumiText += "}";
  }
  else if ( iPeriod==12 ) {
    lumiText += "8 TeV";
  }
  else if ( iPeriod==0 ) {
    lumiText += lumi_sqrtS;
  }
   
  std::cout << lumiText << std::endl;

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(lumiTextSize*t);    
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  if( outOfFrame ) {
    latex.SetTextFont(cmsTextFont);
    latex.SetTextAlign(11); 
    latex.SetTextSize(cmsTextSize*t);    
    latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);
  }
  
  pad->cd();

  float posX_=0;
  if( iPosX%10<=1 ) {
    posX_ =   l + relPosX*(1-l-r);
  }
  else if( iPosX%10==2 ) {
    posX_ =  l + 0.5*(1-l-r);
  }
  else if( iPosX%10==3 ) {
    posX_ =  1-r - relPosX*(1-l-r);
  }
  float posY_ = 1-t - relPosY*(1-t-b);

  if( !outOfFrame ) {
    if( drawLogo ) {
      posX_ =   l + 0.045*(1-l-r)*W/H;
      posY_ = 1-t - 0.045*(1-t-b);
      float xl_0 = posX_;
      float yl_0 = posY_ - 0.15;
      float xl_1 = posX_ + 0.15*H/W;
      float yl_1 = posY_;
      TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
      TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
      pad_logo->Draw();
      pad_logo->cd();
      CMS_logo->Draw("X");
      pad_logo->Modified();
      pad->cd();
    } else {
      latex.SetTextFont(cmsTextFont);
      latex.SetTextSize(cmsTextSize*t);
      latex.SetTextAlign(align_);
      latex.DrawLatex(posX_, posY_, cmsText);
      if( writeExtraText ) {
        latex.SetTextFont(extraTextFont);
        latex.SetTextAlign(align_);
        latex.SetTextSize(extraTextSize*t);
        latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
      }
    }
  }
  else if( writeExtraText ) {
    if( iPosX==0) {
      posX_ =   l +  relPosX*(1-l-r);
      posY_ =   1-t+lumiTextOffset*t;
    }
    latex.SetTextFont(extraTextFont);
    latex.SetTextSize(extraTextSize*t);
    latex.SetTextAlign(align_);
    latex.DrawLatex(posX_, posY_, extraText);      
  }

  return;

}

void LDMX_lumi(TPad* pad, int iPosX, TString extraText)
{

  // Global variables

  TString cmsText     = "LDMX";
  float cmsTextFont   = 52;  // default is helvetic-italics

  bool writeExtraText = true;
  float extraTextFont = 52;  // default is helvetica-italics
  
  // text sizes and text offsets with respect to the top frame
  // in unit of the top margin size
  float lumiTextSize     = 0.6;
  float lumiTextOffset   = 0.2;
  float cmsTextSize      = 0.6;
  float cmsTextOffset    = 0.17; //0.1;  // only used in outOfFrame version

  float relPosX    = 0.3;
  float relPosY    = 0.035;
  float relExtraDY = 1.2;

  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;

  bool drawLogo      = false;

  bool outOfFrame    = false;
  if( iPosX/10==0 ) 
    {
      outOfFrame = true;
    }
  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  if( iPosX == 0  ) relPosX = cmsTextOffset;
  int align_ = 10*alignX_ + alignY_;

  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  //  float e = 0.025;

  pad->cd();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  if( outOfFrame ) {
    latex.SetTextFont(cmsTextFont);
    latex.SetTextAlign(11); 
    latex.SetTextSize(cmsTextSize*t);    
    latex.DrawLatex(l + 0.55*(1-l-r),1-t+lumiTextOffset*t,cmsText);
  }
  
  pad->cd();

  float posX_=0;
  if( iPosX%10<=1 ) {
    posX_ =   l + relPosX*(1-l-r);
  }
  else if( iPosX%10==2 ) {
    posX_ =  l + 0.5*(1-l-r);
  }
  else if( iPosX%10==3 ) {
    posX_ =  1-r - relPosX*(1-l-r);
  }
  float posY_ = 1-t - relPosY*(1-t-b);

  if( !outOfFrame ) {
      latex.SetTextFont(cmsTextFont);
      latex.SetTextSize(cmsTextSize*t);
      latex.SetTextAlign(align_);
      latex.DrawLatex(posX_, posY_, cmsText);
      if( writeExtraText ) {
        latex.SetTextFont(extraTextFont);
        latex.SetTextAlign(align_);
        latex.SetTextSize(extraTextSize*t);
        latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
      }
  }
  else if( writeExtraText ) {
    if( iPosX==0) {
      posX_ =   l +  (0.55+relPosX)*(1-l-r);
      posY_ =   1-t+lumiTextOffset*t;
    }
    latex.SetTextFont(extraTextFont);
    latex.SetTextSize(extraTextSize*t);
    latex.SetTextAlign(align_);
    latex.DrawLatex(posX_, posY_, extraText);      
  }

  return;

}

ColorMap DefaultColors()
{

  ColorMap colors;

  colors["ttbar"]           = color_ttbar;
  colors["ttbar1l"]         = color_ttbar1l;
  colors["ttbar2l"]         = color_ttbar2l;
  colors["ttZ"]             = color_ttZ;
  colors["wjets"]           = color_wjets;
  colors["other"]           = color_other; 
  colors["znunu"]           = color_znunu;
  colors["tW"]              = color_tW;
  colors["ttW"]             = color_ttW;
  colors["qcd"]             = color_qcd;
  colors["T2tt_850_100"]    = color_T2tt_850_100;
  colors["T2tt_650_325"]    = color_T2tt_650_325;
  colors["T2tt_500_325"]    = color_T2tt_500_325;
  colors["T2tt_425_325"]    = color_T2tt_425_325;
  colors["T1tttt_1500_100"] = color_T1tttt_1500_100;
  colors["T1tttt_1200_800"] = color_T1tttt_1200_800;
  colors["comp1"]           = color_comp1;
  colors["comp2"]           = color_comp2;
  colors["comp3"]           = color_comp3;
  colors["comp4"]           = color_comp4;
  colors["comp5"]           = color_comp5;
  colors["comp6"]           = color_comp6;
  colors["comp7"]           = color_comp7;
  colors["comp8"]           = color_comp8;
  colors["comp9"]           = color_comp9;
  colors["comp10"]          = color_comp10;
  colors["comp11"]          = color_comp11;
  colors["comp12"]          = color_comp12;

  return colors;

}

TCanvas* newAlbum(Int_t numItems, const Char_t* name, const Char_t* title, Int_t width, Int_t height, Int_t rows, Int_t columns)
{
  if (rows > 0 && columns > 0) {
  } else if (rows > 0) {
    columns   = TMath::CeilNint(1.0 * numItems / rows);
  } else if (columns > 0) {
    rows      = TMath::CeilNint(1.0 * numItems / columns);
  } else {
    columns   = TMath::CeilNint(TMath::Sqrt(numItems));
    rows      = TMath::CeilNint(1.0 * numItems / columns);
  }
  if (width < 0) {
    if      (columns == 1)  width   = 550;
    else if (columns == 2)  width   = 1100;
    else if (columns == 3)  width   = 1500;
    else                    width   = 1780;
  }
  if (height < 0) {
    if      (rows == 1)     height  = 550;
    else if (rows == 2)     height  = 950;
    else                    height  = 1000;
  }
  TCanvas*    canvas = MakeCanvas(name, title, width, height);
  if (columns > 1 || rows > 1)
    canvas->Divide(columns, rows, 0.00001f, 0.00001f);
  return canvas;
}

TCanvas* drawAll(TObjArray* plots, TString name)
{
  if (plots == 0 || plots->GetEntries() == 0) {
    std::cout << "WARNING : Nothing to draw." << std::endl;
    return 0;
  }

  Int_t             numCols   = 0;
  Int_t             numRows   = 0;


  Int_t             numPlots  = plots->GetEntriesFast();
  TCanvas*          canvas    = newAlbum(numPlots, name, 0, -1, -1, numRows, numCols);
  for (Int_t iPlot = 0; iPlot < numPlots; ++iPlot) {
    TObject*        object    = plots->At(iPlot);
    if (object == 0)          continue;

    TVirtualPad*    pad       = canvas->cd(iPlot + 1);

    if (object->InheritsFrom("TH1")) {
      TH1*        histo     = dynamic_cast<TH1*>(object);
      if (histo->GetMaximum() > histo->GetMinimum()) {
        if (histo->GetDimension() > 1 && pad->GetRightMargin() < 0.1)
          pad->SetRightMargin(0.2f);

        //-- 2D histograms ----------------------------------------------------
        histo->Draw();

        if (pad->GetLogx())
          histo->GetXaxis()->SetMoreLogLabels();
        if (pad->GetLogy() && histo->GetMaximum() < 1e2*histo->GetMinimum())
          histo->GetYaxis()->SetMoreLogLabels();
        if (pad->GetLogz() && histo->GetMaximum() < 1e2*histo->GetMinimum())
          histo->GetZaxis()->SetMoreLogLabels();
      }
    }
    else if (object->InheritsFrom("TCanvas")) {
      TCanvas*      sub     = dynamic_cast<TCanvas*>(object);
      sub->DrawClonePad();
    } else {
      object->Draw();
    }
  }

  canvas->Update();
  return canvas;
}

int colorGetter(int counter){
  int color = counter+1;
  if(color>= 10) color++;
  return color;
}

TH1D* addOverFlow(TH1D* h)
{

  int nbins = h->GetNbinsX()+1;

  double e1 = h->GetBinError(nbins-1);
  double e2 = h->GetBinError(nbins);
  h->AddBinContent(nbins-1, h->GetBinContent(nbins));
  h->SetBinError(nbins-1, sqrt(e1*e1 + e2*e2));
  h->SetBinContent(nbins, 0);
  h->SetBinError(nbins, 0);
  return h;
}

};

#endif

