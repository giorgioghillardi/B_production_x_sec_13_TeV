#include <TCanvas.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TPave.h>
#include <TLine.h>
#include <TStyle.h>
#include <sstream>

TCanvas *canvasDressing(TString name = "c1")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineWidth(2.0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetErrorX(0.5);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetLabelFont(42,"z");
  gStyle->SetTitleFont(42,"x");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetTitleFont(42,"z");
  gStyle->SetEndErrorSize(1);
    
  TCanvas *c1 = new TCanvas(name,name,800,600);
  //c1->Range(4.819967,-13.30288,6.070376,78.37912);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  // c1->SetLeftMargin(0.1439791);
  // c1->SetRightMargin(0.03534031);
  // c1->SetTopMargin(0.07777778);
  // c1->SetBottomMargin(0.145098);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
    
  return c1;
}

void Legend(int channel, int pt_min, int pt_max, double y_min, double y_max, bool pt_bins_flag)
{
  double y_top = 0.92;
  double x_top_1 = 0.1;
  double x_top_2 = 0.7;

  double x_pt_lable = 0.13;
  double y_pt_lable = 0.84;

  TLatex * tex = new TLatex(x_top_2,y_top,"2.71 fb^{-1} (13 TeV)");
  tex->SetNDC(kTRUE);
  tex->SetLineWidth(2);
  tex->Draw();
  tex = new TLatex(x_top_1,y_top,"CMS Preliminary");
  tex->SetNDC(kTRUE);
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->Draw();
  
  if(pt_bins_flag)
    tex = new TLatex(x_pt_lable, y_pt_lable, TString::Format("%d < p_{T} < %d GeV, %.2f < |y| < %.2f",pt_min, pt_max, y_min, y_max));

  tex->SetNDC(kTRUE);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.04);
  tex->Draw();
}
