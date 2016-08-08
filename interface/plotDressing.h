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
//#########################################
void Legend(int channel, int pt_low, int pt_high, bool pt_bins_flag)
{
double x_pos = 0.65;
  double y_pos_s = 0.93;
  double y_pos = 0.91;
  double y_space = -0.075;
  double x1_line = x_pos-0.07;
  double x2_line = x_pos-0.02;
  double y_line = y_pos+0.02;
  double x_signal = x_pos + 0.1;
  double y_separation = -0.05;

    TPaveText *pt = new TPaveText(0.3,0.9,0.6,0.95,"blNDC");
    pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextFont(42);
    pt->Draw();
    TLatex * tex = new TLatex(x_pos,y_pos_s,"#surds = 13 TeV");
    tex->SetNDC(kTRUE);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(0.1,y_pos_s,"CMS Preliminary");
    tex->SetNDC(kTRUE);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    /*tex = new TLatex(x_pos,y_pos + 2*y_space,"Preliminary");
    tex->SetNDC(kTRUE);
    tex->SetTextFont(52);
    tex->SetLineWidth(2);
    tex->Draw();*/

    if(pt_bins_flag)
      tex = new TLatex(0.55,y_pos + y_space-0.04,TString::Format("%d < p_{T} < %d [GeV]",pt_low,pt_high));
    else
      tex = new TLatex(0.55,y_pos + y_space-0.04 ,TString::Format("p_{T} > 10 GeV"));}
//###########################################

 tex->SetNDC(kTRUE);
    tex->SetLineWidth(2);
    tex->Draw();

    TLine *line = new TLine(x1_line,y_separation+y_line+4*y_space,x2_line,y_separation+y_line+4*y_space);
    line->SetNDC(kTRUE);
    line->SetLineColor(4);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(x_pos,y_separation+y_pos+4*y_space,"Total fit");
    tex->SetNDC(kTRUE);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    line = new TLine(x1_line,y_separation+y_line+5*y_space,x2_line,y_separation+y_line+5*y_space);
    line->SetNDC(kTRUE);
    line->SetLineColor(9);
    line->SetLineStyle(2);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(x_pos,y_separation+y_pos+5*y_space,"Background");
    tex->SetNDC(kTRUE);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    TPave *pave = new TPave(x1_line,y_separation+y_line+6.75*y_space,x2_line,y_separation+y_line+7.25*y_space,4,$
    pave->SetFillColor(8);
    pave->SetFillStyle(3008);
    pave->SetLineColor(30);
    pave->Draw();
    tex = new TLatex(x_pos,y_separation+y_pos+7*y_space,"Signal");
    tex->SetNDC(kTRUE);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();

switch(channel){
    case 1:
      tex = new TLatex(x_pos,y_pos + 3*y_space,"B^{#pm} #rightarrow J/#psi K^{#pm}");
      tex->SetNDC(kTRUE);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(x_pos,y_separation+y_pos+2*y_space,"Peaking background");
      tex->SetNDC(kTRUE);
      tex->SetTextFont(42);
      tex->SetLineWidth(2);
      tex->Draw();
      line = new TLine(x1_line,y_separation+y_line+7*y_space,x2_line,y_separation+y_line+7*y_space);
      line->SetNDC(kTRUE); 
      line->SetLineColor(6);
      line->SetLineStyle(2);
      line->SetLineWidth(3);
      line->Draw();
      break;
    case 2:
      tex = new TLatex(0.55,y_pos + 2*y_space-0.04,"B^{0}_{d} #rightarrow J/#psi K^{0*}");
      tex->SetNDC(kTRUE);
      tex->SetLineWidth(2);
      tex->Draw();
      break;
    case 3:
      tex = new TLatex(x_pos,y_pos + 2*y_space,"B^{0}_{d} #rightarrow J/#psi K_{s}");
      tex->SetNDC(kTRUE);
      tex->SetLineWidth(2);
      tex->Draw();
      break;
    case 4:
      tex = new TLatex(0.55,y_pos + 2*y_space-0.04,"B^{0}_{s} #rightarrow J/#psi #phi");
      tex->SetNDC(kTRUE);
      tex->SetLineWidth(2);
      tex->Draw();
      break;
    case 5:
      tex = new TLatex(x_pos,y_pos + 3*y_space,"#psi_{2S}#rightarrowJ/#psi #pi^{+}#pi^{-}");
      tex->SetNDC(kTRUE);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(x_signal,y_separation+y_pos+8*y_space,"(#psi_{2S})");
      tex->SetNDC(kTRUE);
      tex->SetTextFont(42);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(x_pos,y_separation+y_pos+9*y_space,"Signal");
      tex->SetNDC(kTRUE);
      tex->SetTextFont(42);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(x_signal,y_separation+y_pos+9*y_space,"(X(3872))");
      tex->SetNDC(kTRUE);
      tex->SetTextFont(42);
      tex->SetLineWidth(2);
      tex->Draw();
      pave = new TPave(x1_line,y_separation+y_line+8.75*y_space,x2_line,y_separation+y_line+9.25*y_space,4,"blND$
      pave->SetFillColor(kOrange);
      pave->SetFillStyle(3008);
      pave->SetLineColor(2);
      pave->Draw();
      break;
    case 6:
      tex = new TLatex(x_pos,y_pos + 3*y_space,"#Lambda_{b} #rightarrow J/#psi #Lambda_{0}");
      tex->SetNDC(kTRUE);
      tex->SetLineWidth(2);
      tex->Draw();
      break;
    }






void LegendChannelOne()
{
  double x_pos = 0.65;
  double y_pos = 0.93;
  double y_space = -0.075;
  double x1_line = x_pos-0.07;
  double x2_line = x_pos-0.02;
  double y_line = y_pos+0.02;

    TPaveText *pt = new TPaveText(0.3,0.9,0.6,0.95,"blNDC");
    pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextFont(42);
    pt->Draw();
    TLatex * tex = new TLatex(x_pos,y_pos,"#surds = 13 TeV");
    tex->SetNDC(kTRUE);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(x_pos,y_pos + y_space,"CMS");
    tex->SetNDC(kTRUE);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(x_pos,y_pos + 2*y_space,"Preliminary");
    tex->SetNDC(kTRUE);
    tex->SetTextFont(52);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(x_pos,y_pos + 3*y_space,"B^{+} #rightarrow J/#psi K^{+}");
    tex->SetNDC(kTRUE);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(x_pos,y_pos + 4*y_space,"p_{T} > 10 GeV");
    tex->SetNDC(kTRUE);
    tex->SetLineWidth(2);
    tex->Draw();
    TLine *line = new TLine(x1_line,y_line+6*y_space,x2_line,y_line+6*y_space);
    line->SetNDC(kTRUE);
    line->SetLineColor(4);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(x_pos,y_pos + 6*y_space,"Total fit");
    tex->SetNDC(kTRUE);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    line = new TLine(x1_line,y_line+7*y_space,x2_line,y_line+7*y_space);
    line->SetNDC(kTRUE);
    line->SetLineColor(7);
    line->SetLineStyle(2);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(x_pos,y_pos + 7*y_space,"Background");
    tex->SetNDC(kTRUE);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(x_pos,y_pos + 8*y_space,"Peaking background");
    tex->SetNDC(kTRUE);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    line = new TLine(x1_line,y_line+8*y_space,x2_line,y_line+8*y_space);
    line->SetNDC(kTRUE); 
    line->SetLineColor(6);
    line->SetLineStyle(2);
    line->SetLineWidth(3);
    line->Draw();
    TPave *pave = new TPave(x1_line,y_line+8.75*y_space,x2_line,y_line+9.25*y_space,4,"blNDC");
    pave->SetFillColor(2);
    pave->SetFillStyle(3008);
    pave->SetLineColor(2);
    pave->Draw();
    tex = new TLatex(x_pos,y_pos + 9*y_space,"Signal");
    tex->SetNDC(kTRUE);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
}

void LegendChannelTwo()
{
    TPaveText *pt = new TPaveText(0.3347987,0.94,0.6652013,0.995,"blNDC");
    pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextFont(42);
    pt->Draw();
    TLatex *   tex = new TLatex(5.787879,810.9766,"#surds = 13 TeV");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(5.814992,714.2618,"CMS");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(5.703349,650,"Preliminary");
    tex->SetTextFont(52);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(5.706539,590.0,"B^{0} #rightarrow J/#psi K*");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(5.708134,541.1741,"p_{T} > 10 GeV");
    tex->SetLineWidth(2);
    tex->Draw();
    TLine *line = new TLine(5.637959,488.5939,5.679426,488.5939);
    line->SetLineColor(4);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(5.706539,475.3193,"Total fit");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    line = new TLine(5.642743,427.9101,5.684211,427.9101);
    line->SetLineColor(7);
    line->SetLineStyle(2);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(5.711324,420.3246,"Background");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    TPave *pave = new TPave(5.639553,346.3662,5.685805,386.19,4,"br");
    pave->SetFillColor(2);
    pave->SetFillStyle(3008);
    pave->SetLineColor(2);
    pave->Draw();
    tex = new TLatex(5.711324,352.0553,"Signal");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    
}

void LegendChannelFour()
{
    TPaveText *pt = new TPaveText(0.3347987,0.94,0.6652013,0.995,"blNDC");
    pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextFont(42);
    pt->Draw();
    TLatex *   tex = new TLatex(5.794258,217.5708,"#surds = 13 TeV");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(5.818182,191.2995,"CMS");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(5.706539,176.6482,"Preliminary");
    tex->SetTextFont(52);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(5.706539,159.9759,"B^{0}_{s} #rightarrow J/#psi #phi");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(5.706539,142.2933,"p_{T} > 10 GeV");
    tex->SetLineWidth(2);
    tex->Draw();
    TLine *line = new TLine(5.639553,125.6211,5.681021,125.6211);
    line->SetLineColor(4);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(5.706539,123.6002,"Total fit");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    line = new TLine(5.637959,110.4645,5.682616,110.4645);
    line->SetLineColor(7);
    line->SetLineStyle(2);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(5.706539,107.9384,"Background");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    TPave *pave = new TPave(5.639553,92.0,5.685805,100.0,4,"br");
    pave->SetFillColor(2);
    pave->SetFillStyle(3008);
    pave->SetLineColor(2);
    pave->Draw();
    tex = new TLatex(5.706539,92.0,"Signal");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    
}

void LegendChannelFive()
{
    TPaveText *pt = new TPaveText(0.3347987,0.94,0.6652013,0.995,"blNDC");
    pt->SetName("title");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextFont(42);
    pt->Draw();
    TLatex *   tex = new TLatex(3.916427,2983.799,"#surds = 13 TeV");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(3.932376,2615.727,"CMS");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(3.888995,2421.273,"Preliminary");
    tex->SetTextFont(52);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(3.888995,2171.262,"#psi_{2S}#rightarrowJ/#psi #pi^{+}#pi^{-}");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(3.888995,1886.526,"p_{T} > 10 GeV");
    tex->SetLineWidth(2);
    tex->Draw();
    TLine *line = new TLine(3.857735,1622.625,3.873684,1622.625);
    line->SetLineColor(4);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(3.888995,1615.68,"Total fit");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    line = new TLine(3.859649,1497.619,3.875598,1497.619);
    line->SetLineColor(7);
    line->SetLineStyle(2);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(3.888995,1394.576,"Background");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    TPave *pave = new TPave(3.861659,462.0314,3.876056,547.3295,4,"br");
    pave->SetFillColor(2);
    pave->SetFillStyle(3008);
    pave->SetLineColor(2);
    pave->Draw();
    tex = new TLatex(3.888995,483.3559,"Signal_{#psi_{2S}}");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(3.903599,263.0025,"Signal(X(3872))");
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    pave = new TPave(3.864163,263.0025,3.877934,348.3006,4,"br");
    pave->SetFillColor(kOrange);
    pave->SetFillStyle(3008);
    pave->SetLineColor(2);
    pave->Draw();
    
}


