{
//=========Macro generated from canvas: c1/c1
//=========  (Sun May 11 01:38:19 2014) by ROOT version5.32/00
   TCanvas *c1 = new TCanvas("c1", "c1",294,80,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->SetHighLightColor(2);
   c1->Range(-7.297297,-0.3355263,33.24324,1.638158);
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
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderSize(0);
   
   TGraphErrors *gre = new TGraphErrors(6);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#660066");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#660066");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(20);
   gre->SetPoint(0,2.5,0.8108289);
   gre->SetPointError(0,0,0.02328085);
   gre->SetPoint(1,7.5,0.7438325);
   gre->SetPointError(1,0,0.01032901);
   gre->SetPoint(2,12.5,0.6567098);
   gre->SetPointError(2,0,0.006426298);
   gre->SetPoint(3,17.5,0.5466152);
   gre->SetPointError(3,0,0.006192375);
   gre->SetPoint(4,22.5,0.4512079);
   gre->SetPointError(4,0,0.009264444);
   gre->SetPoint(5,27.5,0.3660221);
   gre->SetPointError(5,0,0.02248456);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,0,30);
   Graph_Graph1->SetMinimum(0);
   Graph_Graph1->SetMaximum(1.5);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->GetXaxis()->SetTitle("N_{Vtx}");
   Graph_Graph1->GetXaxis()->SetNdivisions(505);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetTitle("jet veto efficiency");
   Graph_Graph1->GetYaxis()->SetNdivisions(505);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleOffset(1.5);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetNdivisions(707);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1);
   
   gre->Draw("alp");
   
   gre = new TGraphErrors(6);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#3300cc");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#3300cc");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(21);
   gre->SetPoint(0,2.5,0.9304813);
   gre->SetPointError(0,0,0.02493952);
   gre->SetPoint(1,7.5,0.9323006);
   gre->SetPointError(1,0,0.01156377);
   gre->SetPoint(2,12.5,0.9230914);
   gre->SetPointError(2,0,0.007618974);
   gre->SetPoint(3,17.5,0.909786);
   gre->SetPointError(3,0,0.007988882);
   gre->SetPoint(4,22.5,0.8877687);
   gre->SetPointError(4,0,0.01299514);
   gre->SetPoint(5,27.5,0.8618785);
   gre->SetPointError(5,0,0.03450275);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,0,30);
   Graph_Graph2->SetMinimum(0.8145712);
   Graph_Graph2->SetMaximum(0.9682253);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);
   Graph_Graph2->GetXaxis()->SetNdivisions(505);
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph2->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetNdivisions(505);
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph2->GetYaxis()->SetTitleOffset(1.5);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetNdivisions(707);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph2->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph2->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph2);
   
   gre->Draw("lp");
   
   gre = new TGraphErrors(6);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#660066");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#660066");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(24);
   gre->SetPoint(0,2.5,0.8402406);
   gre->SetPointError(0,0,0.02369933);
   gre->SetPoint(1,7.5,0.8432301);
   gre->SetPointError(1,0,0.01099751);
   gre->SetPoint(2,12.5,0.8332285);
   gre->SetPointError(2,0,0.007238627);
   gre->SetPoint(3,17.5,0.8183795);
   gre->SetPointError(3,0,0.007576939);
   gre->SetPoint(4,22.5,0.7998859);
   gre->SetPointError(4,0,0.01233517);
   gre->SetPoint(5,27.5,0.7679558);
   gre->SetPointError(5,0,0.03256858);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","Graph",100,0,30);
   Graph_Graph3->SetMinimum(0.7225319);
   Graph_Graph3->SetMaximum(0.8767952);
   Graph_Graph3->SetDirectory(0);
   Graph_Graph3->SetStats(0);
   Graph_Graph3->GetXaxis()->SetNdivisions(505);
   Graph_Graph3->GetXaxis()->SetLabelFont(42);
   Graph_Graph3->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph3->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph3->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph3->GetXaxis()->SetTitleFont(42);
   Graph_Graph3->GetYaxis()->SetNdivisions(505);
   Graph_Graph3->GetYaxis()->SetLabelFont(42);
   Graph_Graph3->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph3->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph3->GetYaxis()->SetTitleOffset(1.5);
   Graph_Graph3->GetYaxis()->SetTitleFont(42);
   Graph_Graph3->GetZaxis()->SetNdivisions(707);
   Graph_Graph3->GetZaxis()->SetLabelFont(42);
   Graph_Graph3->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph3->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph3->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph3->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph3->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph3);
   
   gre->Draw("lp");
   
   gre = new TGraphErrors(6);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#3300cc");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#3300cc");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(25);
   gre->SetPoint(0,2.5,0.9385027);
   gre->SetPointError(0,0,0.02504678);
   gre->SetPoint(1,7.5,0.9434882);
   gre->SetPointError(1,0,0.01163294);
   gre->SetPoint(2,12.5,0.9460445);
   gre->SetPointError(2,0,0.007713117);
   gre->SetPoint(3,17.5,0.9461943);
   gre->SetPointError(3,0,0.008147166);
   gre->SetPoint(4,22.5,0.9368461);
   gre->SetPointError(4,0,0.0133495);
   gre->SetPoint(5,27.5,0.9350829);
   gre->SetPointError(5,0,0.03593815);
   
   TH1F *Graph_Graph4 = new TH1F("Graph_Graph4","Graph",100,0,30);
   Graph_Graph4->SetMinimum(0.8919571);
   Graph_Graph4->SetMaximum(0.9782087);
   Graph_Graph4->SetDirectory(0);
   Graph_Graph4->SetStats(0);
   Graph_Graph4->GetXaxis()->SetNdivisions(505);
   Graph_Graph4->GetXaxis()->SetLabelFont(42);
   Graph_Graph4->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph4->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph4->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph4->GetXaxis()->SetTitleFont(42);
   Graph_Graph4->GetYaxis()->SetNdivisions(505);
   Graph_Graph4->GetYaxis()->SetLabelFont(42);
   Graph_Graph4->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph4->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph4->GetYaxis()->SetTitleOffset(1.5);
   Graph_Graph4->GetYaxis()->SetTitleFont(42);
   Graph_Graph4->GetZaxis()->SetNdivisions(707);
   Graph_Graph4->GetZaxis()->SetLabelFont(42);
   Graph_Graph4->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph4->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph4->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph4->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph4->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph4);
   
   gre->Draw("lp");
   
   TLegend *leg = new TLegend(0.55,0.77,0.8,0.85,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("Graph","p_{T} > 20 GeV","lp");

   ci = TColor::GetColor("#660066");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#660066");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph","p_{T} > 30 GeV","lp");

   ci = TColor::GetColor("#3300cc");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#3300cc");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   TLatex *   tex = new TLatex(0.2,0.87,"ATLAS");
tex->SetNDC();
   tex->SetTextFont(72);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.33,0.87,"Simulation Preliminary");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.82,"qq' #rightarrow Hqq', H#rightarrowZZ#rightarrow 4l");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.78,"Anti-k_{t} LCW+JES R=0.4");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.74,"solid markers: inclusive");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.7,"open markers: JVT>0.2");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
