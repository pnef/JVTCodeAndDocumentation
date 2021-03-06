{
//=========Macro generated from canvas: c1/Histogram Drawing Options
//=========  (Sun May 11 01:27:49 2014) by ROOT version5.32/00
   TCanvas *c1 = new TCanvas("c1", "Histogram Drawing Options",193,139,600,800);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->SetHighLightColor(2);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.18);
   c1->SetRightMargin(0.08);
   c1->SetTopMargin(0.07);
   c1->SetBottomMargin(0.17);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderSize(0);
  
// ------------>Primitives in pad: pad1
   TPad *pad1 = new TPad("pad1", "",0,0.25,1,1);
   pad1->Draw();
   pad1->cd();
   pad1->Range(-3.567568,0.7328947,2.918919,1.127632);
   pad1->SetFillColor(0);
   pad1->SetBorderMode(0);
   pad1->SetBorderSize(1);
   pad1->SetTickx(1);
   pad1->SetTicky(1);
   pad1->SetLeftMargin(0.18);
   pad1->SetRightMargin(0.08);
   pad1->SetTopMargin(0.07);
   pad1->SetBottomMargin(0.17);
   pad1->SetFrameBorderMode(0);
   pad1->SetFrameBorderSize(0);
   pad1->SetFrameBorderMode(0);
   pad1->SetFrameBorderSize(0);
   
   TH1D *EffRecoCorrSystSmooth = new TH1D("EffRecoCorrSystSmooth","",10,-2.4,2.4);
   EffRecoCorrSystSmooth->SetBinContent(1,0.9042811);
   EffRecoCorrSystSmooth->SetBinContent(2,0.9041257);
   EffRecoCorrSystSmooth->SetBinContent(3,0.9095281);
   EffRecoCorrSystSmooth->SetBinContent(4,0.9268129);
   EffRecoCorrSystSmooth->SetBinContent(5,0.9236186);
   EffRecoCorrSystSmooth->SetBinContent(6,0.9100666);
   EffRecoCorrSystSmooth->SetBinContent(7,0.927394);
   EffRecoCorrSystSmooth->SetBinContent(8,0.9130248);
   EffRecoCorrSystSmooth->SetBinContent(9,0.9047548);
   EffRecoCorrSystSmooth->SetBinContent(10,0.9048366);
   EffRecoCorrSystSmooth->SetBinError(1,0.01280908);
   EffRecoCorrSystSmooth->SetBinError(2,0.01614723);
   EffRecoCorrSystSmooth->SetBinError(3,0.01162822);
   EffRecoCorrSystSmooth->SetBinError(4,0.006571851);
   EffRecoCorrSystSmooth->SetBinError(5,0.01190353);
   EffRecoCorrSystSmooth->SetBinError(6,0.01194951);
   EffRecoCorrSystSmooth->SetBinError(7,0.006902261);
   EffRecoCorrSystSmooth->SetBinError(8,0.01166346);
   EffRecoCorrSystSmooth->SetBinError(9,0.01628981);
   EffRecoCorrSystSmooth->SetBinError(10,0.01341477);
   EffRecoCorrSystSmooth->SetMinimum(0.8);
   EffRecoCorrSystSmooth->SetMaximum(1.1);
   EffRecoCorrSystSmooth->SetEntries(32005);
   EffRecoCorrSystSmooth->SetStats(0);
   EffRecoCorrSystSmooth->SetFillColor(17);
   EffRecoCorrSystSmooth->SetLineColor(17);
   EffRecoCorrSystSmooth->SetMarkerColor(17);
   EffRecoCorrSystSmooth->GetXaxis()->SetTitle("jet #eta");
   EffRecoCorrSystSmooth->GetXaxis()->SetNdivisions(505);
   EffRecoCorrSystSmooth->GetXaxis()->SetLabelFont(42);
   EffRecoCorrSystSmooth->GetXaxis()->SetLabelSize(0.05);
   EffRecoCorrSystSmooth->GetXaxis()->SetTitleSize(0.05);
   EffRecoCorrSystSmooth->GetXaxis()->SetTitleOffset(1.1);
   EffRecoCorrSystSmooth->GetXaxis()->SetTitleFont(42);
   EffRecoCorrSystSmooth->GetYaxis()->SetTitle("Efficiency");
   EffRecoCorrSystSmooth->GetYaxis()->SetNdivisions(505);
   EffRecoCorrSystSmooth->GetYaxis()->SetLabelFont(42);
   EffRecoCorrSystSmooth->GetYaxis()->SetLabelSize(0.05);
   EffRecoCorrSystSmooth->GetYaxis()->SetTitleSize(0.05);
   EffRecoCorrSystSmooth->GetYaxis()->SetTitleOffset(1.5);
   EffRecoCorrSystSmooth->GetYaxis()->SetTitleFont(42);
   EffRecoCorrSystSmooth->GetZaxis()->SetNdivisions(707);
   EffRecoCorrSystSmooth->GetZaxis()->SetLabelFont(42);
   EffRecoCorrSystSmooth->GetZaxis()->SetLabelOffset(0.015);
   EffRecoCorrSystSmooth->GetZaxis()->SetLabelSize(0.05);
   EffRecoCorrSystSmooth->GetZaxis()->SetTitleSize(0.05);
   EffRecoCorrSystSmooth->GetZaxis()->SetTitleOffset(1.1);
   EffRecoCorrSystSmooth->GetZaxis()->SetTitleFont(42);
   EffRecoCorrSystSmooth->Draw("E2");
   
   TH1D *EffRecoCorr = new TH1D("EffRecoCorr","",10,-2.4,2.4);
   EffRecoCorr->SetBinContent(1,0.9042811);
   EffRecoCorr->SetBinContent(2,0.9041257);
   EffRecoCorr->SetBinContent(3,0.9095281);
   EffRecoCorr->SetBinContent(4,0.9268129);
   EffRecoCorr->SetBinContent(5,0.9236186);
   EffRecoCorr->SetBinContent(6,0.9100666);
   EffRecoCorr->SetBinContent(7,0.927394);
   EffRecoCorr->SetBinContent(8,0.9130248);
   EffRecoCorr->SetBinContent(9,0.9047548);
   EffRecoCorr->SetBinContent(10,0.9048366);
   EffRecoCorr->SetBinError(1,0.008415822);
   EffRecoCorr->SetBinError(2,0.007999708);
   EffRecoCorr->SetBinError(3,0.007135206);
   EffRecoCorr->SetBinError(4,0.005873792);
   EffRecoCorr->SetBinError(5,0.006147647);
   EffRecoCorr->SetBinError(6,0.006474346);
   EffRecoCorr->SetBinError(7,0.006240395);
   EffRecoCorr->SetBinError(8,0.007147195);
   EffRecoCorr->SetBinError(9,0.00826719);
   EffRecoCorr->SetBinError(10,0.009305642);
   EffRecoCorr->SetEntries(32005);
   EffRecoCorr->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#3300cc");
   EffRecoCorr->SetLineColor(ci);

   ci = TColor::GetColor("#3300cc");
   EffRecoCorr->SetMarkerColor(ci);
   EffRecoCorr->SetMarkerStyle(20);
   EffRecoCorr->GetXaxis()->SetTitle("jet #eta");
   EffRecoCorr->GetXaxis()->SetNdivisions(505);
   EffRecoCorr->GetXaxis()->SetLabelFont(42);
   EffRecoCorr->GetXaxis()->SetLabelSize(0.05);
   EffRecoCorr->GetXaxis()->SetTitleSize(0.05);
   EffRecoCorr->GetXaxis()->SetTitleOffset(1.1);
   EffRecoCorr->GetXaxis()->SetTitleFont(42);
   EffRecoCorr->GetYaxis()->SetTitle("Efficiency");
   EffRecoCorr->GetYaxis()->SetNdivisions(505);
   EffRecoCorr->GetYaxis()->SetLabelFont(42);
   EffRecoCorr->GetYaxis()->SetLabelSize(0.05);
   EffRecoCorr->GetYaxis()->SetTitleSize(0.05);
   EffRecoCorr->GetYaxis()->SetTitleOffset(1.5);
   EffRecoCorr->GetYaxis()->SetTitleFont(42);
   EffRecoCorr->GetZaxis()->SetNdivisions(707);
   EffRecoCorr->GetZaxis()->SetLabelFont(42);
   EffRecoCorr->GetZaxis()->SetLabelOffset(0.015);
   EffRecoCorr->GetZaxis()->SetLabelSize(0.05);
   EffRecoCorr->GetZaxis()->SetTitleSize(0.05);
   EffRecoCorr->GetZaxis()->SetTitleOffset(1.1);
   EffRecoCorr->GetZaxis()->SetTitleFont(42);
   EffRecoCorr->Draw("EXsame");
   
   TH1D *EffDataCorr = new TH1D("EffDataCorr","",10,-2.4,2.4);
   EffDataCorr->SetBinContent(1,0.8963016);
   EffDataCorr->SetBinContent(2,0.9000123);
   EffDataCorr->SetBinContent(3,0.9106385);
   EffDataCorr->SetBinContent(4,0.9202884);
   EffDataCorr->SetBinContent(5,0.9173704);
   EffDataCorr->SetBinContent(6,0.9236579);
   EffDataCorr->SetBinContent(7,0.9234973);
   EffDataCorr->SetBinContent(8,0.9043903);
   EffDataCorr->SetBinContent(9,0.9058377);
   EffDataCorr->SetBinContent(10,0.8930726);
   EffDataCorr->SetBinError(1,0.004389952);
   EffDataCorr->SetBinError(2,0.003782319);
   EffDataCorr->SetBinError(3,0.003253793);
   EffDataCorr->SetBinError(4,0.002922135);
   EffDataCorr->SetBinError(5,0.002937215);
   EffDataCorr->SetBinError(6,0.002870226);
   EffDataCorr->SetBinError(7,0.002937959);
   EffDataCorr->SetBinError(8,0.003408607);
   EffDataCorr->SetBinError(9,0.003677489);
   EffDataCorr->SetBinError(10,0.004472196);
   EffDataCorr->SetEntries(77372);
   EffDataCorr->SetStats(0);
   EffDataCorr->SetMarkerStyle(20);
   EffDataCorr->GetXaxis()->SetNdivisions(505);
   EffDataCorr->GetXaxis()->SetLabelFont(42);
   EffDataCorr->GetXaxis()->SetLabelSize(0.05);
   EffDataCorr->GetXaxis()->SetTitleSize(0.05);
   EffDataCorr->GetXaxis()->SetTitleOffset(1.1);
   EffDataCorr->GetXaxis()->SetTitleFont(42);
   EffDataCorr->GetYaxis()->SetNdivisions(505);
   EffDataCorr->GetYaxis()->SetLabelFont(42);
   EffDataCorr->GetYaxis()->SetLabelSize(0.05);
   EffDataCorr->GetYaxis()->SetTitleSize(0.05);
   EffDataCorr->GetYaxis()->SetTitleOffset(1.5);
   EffDataCorr->GetYaxis()->SetTitleFont(42);
   EffDataCorr->GetZaxis()->SetNdivisions(707);
   EffDataCorr->GetZaxis()->SetLabelFont(42);
   EffDataCorr->GetZaxis()->SetLabelOffset(0.015);
   EffDataCorr->GetZaxis()->SetLabelSize(0.05);
   EffDataCorr->GetZaxis()->SetTitleSize(0.05);
   EffDataCorr->GetZaxis()->SetTitleOffset(1.1);
   EffDataCorr->GetZaxis()->SetTitleFont(42);
   EffDataCorr->Draw("EXsame");
   
   TLegend *leg = new TLegend(0.55,0.73,0.8,0.85,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("EffRecoCorrSystSmooth","Total Uncertainty","fl");
   entry->SetFillColor(17);
   entry->SetFillStyle(1001);
   entry->SetLineColor(17);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("EffRecoCorr","MC","lp");

   ci = TColor::GetColor("#3300cc");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#3300cc");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("EffDataCorr","Data","lp");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   leg->Draw();
   TLatex *   tex = new TLatex(0.2,0.87,"ATLAS");
tex->SetNDC();
   tex->SetTextFont(72);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.33,0.87,"Preliminary");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.82,"#sqrt{s} = 8 TeV, L = 20.3 fb^{-1}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.78,"Sherpa Z#rightarrow#mu#mu");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.74,"Anti-k_{t} LCW+JES R=0.4");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.7,"20 < p_{T} < 50 GeV, |#eta| < 2.4");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.66,"40 < p_{T}^{ref} < 50 GeV");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.62,"JVT > 0.7");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   pad1->Modified();
   c1->cd();
  
// ------------>Primitives in pad: pad2
   pad2 = new TPad("pad2", "",0,0,1,0.25);
   pad2->Draw();
   pad2->cd();
   pad2->Range(-3.567568,0.8996226,2.918919,1.050566);
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(2);
   pad2->SetTickx(1);
   pad2->SetTicky(1);
   pad2->SetLeftMargin(0.18);
   pad2->SetRightMargin(0.08);
   pad2->SetTopMargin(0.07);
   pad2->SetBottomMargin(0.4);
   pad2->SetFrameBorderMode(0);
   pad2->SetFrameBorderSize(0);
   pad2->SetFrameBorderMode(0);
   pad2->SetFrameBorderSize(0);
   
   TH1D *EffRecoCorrSystSmoothuncert__1 = new TH1D("EffRecoCorrSystSmoothuncert__1","",10,-2.4,2.4);
   EffRecoCorrSystSmoothuncert__1->SetBinContent(1,1);
   EffRecoCorrSystSmoothuncert__1->SetBinContent(2,1);
   EffRecoCorrSystSmoothuncert__1->SetBinContent(3,1);
   EffRecoCorrSystSmoothuncert__1->SetBinContent(4,1);
   EffRecoCorrSystSmoothuncert__1->SetBinContent(5,1);
   EffRecoCorrSystSmoothuncert__1->SetBinContent(6,1);
   EffRecoCorrSystSmoothuncert__1->SetBinContent(7,1);
   EffRecoCorrSystSmoothuncert__1->SetBinContent(8,1);
   EffRecoCorrSystSmoothuncert__1->SetBinContent(9,1);
   EffRecoCorrSystSmoothuncert__1->SetBinContent(10,1);
   EffRecoCorrSystSmoothuncert__1->SetBinError(1,0.01416494);
   EffRecoCorrSystSmoothuncert__1->SetBinError(2,0.0178595);
   EffRecoCorrSystSmoothuncert__1->SetBinError(3,0.0127849);
   EffRecoCorrSystSmoothuncert__1->SetBinError(4,0.007090806);
   EffRecoCorrSystSmoothuncert__1->SetBinError(5,0.01288792);
   EffRecoCorrSystSmoothuncert__1->SetBinError(6,0.01313037);
   EffRecoCorrSystSmoothuncert__1->SetBinError(7,0.007442642);
   EffRecoCorrSystSmoothuncert__1->SetBinError(8,0.01277453);
   EffRecoCorrSystSmoothuncert__1->SetBinError(9,0.01800467);
   EffRecoCorrSystSmoothuncert__1->SetBinError(10,0.01482563);
   EffRecoCorrSystSmoothuncert__1->SetMinimum(0.96);
   EffRecoCorrSystSmoothuncert__1->SetMaximum(1.04);
   EffRecoCorrSystSmoothuncert__1->SetEntries(32015);
   EffRecoCorrSystSmoothuncert__1->SetDirectory(0);
   EffRecoCorrSystSmoothuncert__1->SetStats(0);
   EffRecoCorrSystSmoothuncert__1->SetFillColor(15);
   EffRecoCorrSystSmoothuncert__1->SetFillStyle(3001);
   EffRecoCorrSystSmoothuncert__1->SetLineColor(15);
   EffRecoCorrSystSmoothuncert__1->SetMarkerColor(15);
   EffRecoCorrSystSmoothuncert__1->GetXaxis()->SetNdivisions(505);
   EffRecoCorrSystSmoothuncert__1->GetXaxis()->SetLabelFont(42);
   EffRecoCorrSystSmoothuncert__1->GetXaxis()->SetLabelSize(0.15);
   EffRecoCorrSystSmoothuncert__1->GetXaxis()->SetTitleSize(0.15);
   EffRecoCorrSystSmoothuncert__1->GetXaxis()->SetTickLength(0.09);
   EffRecoCorrSystSmoothuncert__1->GetXaxis()->SetTitleOffset(1.1);
   EffRecoCorrSystSmoothuncert__1->GetXaxis()->SetTitleFont(42);
   EffRecoCorrSystSmoothuncert__1->GetYaxis()->SetTitle("Data / MC");
   EffRecoCorrSystSmoothuncert__1->GetYaxis()->SetNdivisions(505);
   EffRecoCorrSystSmoothuncert__1->GetYaxis()->SetLabelFont(42);
   EffRecoCorrSystSmoothuncert__1->GetYaxis()->SetLabelSize(0.15);
   EffRecoCorrSystSmoothuncert__1->GetYaxis()->SetTitleSize(0.15);
   EffRecoCorrSystSmoothuncert__1->GetYaxis()->SetTickLength(0.04);
   EffRecoCorrSystSmoothuncert__1->GetYaxis()->SetTitleOffset(0.51);
   EffRecoCorrSystSmoothuncert__1->GetYaxis()->SetTitleFont(42);
   EffRecoCorrSystSmoothuncert__1->GetZaxis()->SetNdivisions(707);
   EffRecoCorrSystSmoothuncert__1->GetZaxis()->SetLabelFont(42);
   EffRecoCorrSystSmoothuncert__1->GetZaxis()->SetLabelOffset(0.015);
   EffRecoCorrSystSmoothuncert__1->GetZaxis()->SetLabelSize(0.15);
   EffRecoCorrSystSmoothuncert__1->GetZaxis()->SetTitleSize(0.15);
   EffRecoCorrSystSmoothuncert__1->GetZaxis()->SetTitleOffset(1.1);
   EffRecoCorrSystSmoothuncert__1->GetZaxis()->SetTitleFont(42);
   EffRecoCorrSystSmoothuncert__1->Draw("E2");
   TLine *line = new TLine(-2.4,1,2.4,1);
   line->SetLineStyle(7);
   line->SetLineWidth(2);
   line->Draw();
   
   TH1D *EffRecoCorrSystSmoothuncert = new TH1D("EffRecoCorrSystSmoothuncert","",10,-2.4,2.4);
   EffRecoCorrSystSmoothuncert->SetBinContent(1,1);
   EffRecoCorrSystSmoothuncert->SetBinContent(2,1);
   EffRecoCorrSystSmoothuncert->SetBinContent(3,1);
   EffRecoCorrSystSmoothuncert->SetBinContent(4,1);
   EffRecoCorrSystSmoothuncert->SetBinContent(5,1);
   EffRecoCorrSystSmoothuncert->SetBinContent(6,1);
   EffRecoCorrSystSmoothuncert->SetBinContent(7,1);
   EffRecoCorrSystSmoothuncert->SetBinContent(8,1);
   EffRecoCorrSystSmoothuncert->SetBinContent(9,1);
   EffRecoCorrSystSmoothuncert->SetBinContent(10,1);
   EffRecoCorrSystSmoothuncert->SetBinError(1,0.01416494);
   EffRecoCorrSystSmoothuncert->SetBinError(2,0.0178595);
   EffRecoCorrSystSmoothuncert->SetBinError(3,0.0127849);
   EffRecoCorrSystSmoothuncert->SetBinError(4,0.007090806);
   EffRecoCorrSystSmoothuncert->SetBinError(5,0.01288792);
   EffRecoCorrSystSmoothuncert->SetBinError(6,0.01313037);
   EffRecoCorrSystSmoothuncert->SetBinError(7,0.007442642);
   EffRecoCorrSystSmoothuncert->SetBinError(8,0.01277453);
   EffRecoCorrSystSmoothuncert->SetBinError(9,0.01800467);
   EffRecoCorrSystSmoothuncert->SetBinError(10,0.01482563);
   EffRecoCorrSystSmoothuncert->SetMinimum(0.96);
   EffRecoCorrSystSmoothuncert->SetMaximum(1.04);
   EffRecoCorrSystSmoothuncert->SetEntries(32015);
   EffRecoCorrSystSmoothuncert->SetStats(0);
   EffRecoCorrSystSmoothuncert->SetFillColor(15);
   EffRecoCorrSystSmoothuncert->SetFillStyle(3001);
   EffRecoCorrSystSmoothuncert->SetLineColor(15);
   EffRecoCorrSystSmoothuncert->SetMarkerColor(15);
   EffRecoCorrSystSmoothuncert->GetXaxis()->SetNdivisions(505);
   EffRecoCorrSystSmoothuncert->GetXaxis()->SetLabelFont(42);
   EffRecoCorrSystSmoothuncert->GetXaxis()->SetLabelSize(0.15);
   EffRecoCorrSystSmoothuncert->GetXaxis()->SetTitleSize(0.15);
   EffRecoCorrSystSmoothuncert->GetXaxis()->SetTickLength(0.09);
   EffRecoCorrSystSmoothuncert->GetXaxis()->SetTitleOffset(1.1);
   EffRecoCorrSystSmoothuncert->GetXaxis()->SetTitleFont(42);
   EffRecoCorrSystSmoothuncert->GetYaxis()->SetTitle("Data / MC");
   EffRecoCorrSystSmoothuncert->GetYaxis()->SetNdivisions(505);
   EffRecoCorrSystSmoothuncert->GetYaxis()->SetLabelFont(42);
   EffRecoCorrSystSmoothuncert->GetYaxis()->SetLabelSize(0.15);
   EffRecoCorrSystSmoothuncert->GetYaxis()->SetTitleSize(0.15);
   EffRecoCorrSystSmoothuncert->GetYaxis()->SetTickLength(0.04);
   EffRecoCorrSystSmoothuncert->GetYaxis()->SetTitleOffset(0.51);
   EffRecoCorrSystSmoothuncert->GetYaxis()->SetTitleFont(42);
   EffRecoCorrSystSmoothuncert->GetZaxis()->SetNdivisions(707);
   EffRecoCorrSystSmoothuncert->GetZaxis()->SetLabelFont(42);
   EffRecoCorrSystSmoothuncert->GetZaxis()->SetLabelOffset(0.015);
   EffRecoCorrSystSmoothuncert->GetZaxis()->SetLabelSize(0.15);
   EffRecoCorrSystSmoothuncert->GetZaxis()->SetTitleSize(0.15);
   EffRecoCorrSystSmoothuncert->GetZaxis()->SetTitleOffset(1.1);
   EffRecoCorrSystSmoothuncert->GetZaxis()->SetTitleFont(42);
   EffRecoCorrSystSmoothuncert->Draw("E2same");
   
   TH1D *MCsumRatio = new TH1D("MCsumRatio","",10,-2.4,2.4);
   MCsumRatio->SetBinContent(1,0.9911759);
   MCsumRatio->SetBinContent(2,0.9954503);
   MCsumRatio->SetBinContent(3,1.001221);
   MCsumRatio->SetBinContent(4,0.9929603);
   MCsumRatio->SetBinContent(5,0.9932351);
   MCsumRatio->SetBinContent(6,1.014934);
   MCsumRatio->SetBinContent(7,0.9957983);
   MCsumRatio->SetBinContent(8,0.990543);
   MCsumRatio->SetBinContent(9,1.001197);
   MCsumRatio->SetBinContent(10,0.9869988);
   MCsumRatio->SetBinError(1,0.01485555);
   MCsumRatio->SetBinError(2,0.01826381);
   MCsumRatio->SetBinError(3,0.01329101);
   MCsumRatio->SetBinError(4,0.007714584);
   MCsumRatio->SetBinError(5,0.01318985);
   MCsumRatio->SetBinError(6,0.01369458);
   MCsumRatio->SetBinError(7,0.008060053);
   MCsumRatio->SetBinError(8,0.01319296);
   MCsumRatio->SetBinError(9,0.01847879);
   MCsumRatio->SetBinError(10,0.01544506);
   MCsumRatio->SetEntries(50366.93);
   MCsumRatio->SetStats(0);
   MCsumRatio->SetFillColor(1);
   MCsumRatio->SetFillStyle(0);
   MCsumRatio->SetMarkerStyle(20);
   MCsumRatio->GetXaxis()->SetNdivisions(505);
   MCsumRatio->GetXaxis()->SetLabelFont(42);
   MCsumRatio->GetXaxis()->SetLabelSize(0.15);
   MCsumRatio->GetXaxis()->SetTitleSize(0.15);
   MCsumRatio->GetXaxis()->SetTickLength(0.09);
   MCsumRatio->GetXaxis()->SetTitleOffset(1.1);
   MCsumRatio->GetXaxis()->SetTitleFont(42);
   MCsumRatio->GetYaxis()->SetTitle("Data / MC");
   MCsumRatio->GetYaxis()->SetNdivisions(505);
   MCsumRatio->GetYaxis()->SetLabelFont(42);
   MCsumRatio->GetYaxis()->SetLabelSize(0.15);
   MCsumRatio->GetYaxis()->SetTitleSize(0.15);
   MCsumRatio->GetYaxis()->SetTickLength(0.04);
   MCsumRatio->GetYaxis()->SetTitleOffset(0.51);
   MCsumRatio->GetYaxis()->SetTitleFont(42);
   MCsumRatio->GetZaxis()->SetNdivisions(707);
   MCsumRatio->GetZaxis()->SetLabelFont(42);
   MCsumRatio->GetZaxis()->SetLabelOffset(0.015);
   MCsumRatio->GetZaxis()->SetLabelSize(0.15);
   MCsumRatio->GetZaxis()->SetTitleSize(0.15);
   MCsumRatio->GetZaxis()->SetTitleOffset(1.1);
   MCsumRatio->GetZaxis()->SetTitleFont(42);
   MCsumRatio->Draw("EXsame");
   pad2->Modified();
   c1->cd();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
