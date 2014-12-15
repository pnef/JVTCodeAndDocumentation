void SetStyle_MT2Ratio(){ 

// MT2Style ratio ------------------------------------------------------------------------------
  TStyle *MT2StyleRatio= new TStyle("MT2-ratio-Style","MT2-ratio-Style");
  
//  #ifdef __CINT__
//  TStyle *GloStyle;
//  GloStyle = gStyle;                          // save the global style reference
//  gStyle = MT2StyleRatio;
//  #endif
  
  // Canvas
  MT2StyleRatio->SetCanvasColor     (0);
  MT2StyleRatio->SetCanvasBorderMode(0);
  MT2StyleRatio->SetCanvasBorderMode(0);
  MT2StyleRatio->SetCanvasDefH      (200);
  MT2StyleRatio->SetCanvasDefW      (600);

  // Pads
  MT2StyleRatio->SetPadColor       (0);
  MT2StyleRatio->SetPadBorderMode  (0);
  MT2StyleRatio->SetPadBottomMargin(0.4);
  MT2StyleRatio->SetPadTopMargin   (0.07);
  MT2StyleRatio->SetPadLeftMargin  (0.18);
  MT2StyleRatio->SetPadRightMargin (0.08);
  MT2StyleRatio->SetPadTickX          (1);
  MT2StyleRatio->SetPadTickY          (1);

  // Frames
  MT2StyleRatio->SetFrameLineWidth (1);
  MT2StyleRatio->SetFrameFillColor (0);
  MT2StyleRatio->SetFrameBorderMode(0);
  MT2StyleRatio->SetFrameBorderSize(0);


  // Histograms
  MT2StyleRatio->SetHistFillColor(0);
  MT2StyleRatio->SetHistLineWidth(1);

  // Functions
  //MT2StyleRatio->SetFuncColor(1);
  //MT2StyleRatio->SetFuncStyle(0);
  //MT2StyleRatio->SetFuncWidth(2);

  //Legends 
  MT2StyleRatio->SetLegendBorderSize(0);
  MT2StyleRatio->SetFillStyle(0);
  MT2StyleRatio->SetTextFont(42); 

  // Labels, Ticks, and Titles
  MT2StyleRatio->SetTitleSize  ( 0.150,"X");
  MT2StyleRatio->SetTitleOffset( 1.100,"X");
  MT2StyleRatio->SetLabelSize  ( 0.150,"X");
  MT2StyleRatio->SetLabelFont  ( 42   ,"X");
  MT2StyleRatio->SetTitleFont  ( 42   ,"X");
  MT2StyleRatio->SetNdivisions ( 505  ,"X");
  MT2StyleRatio->SetTickLength ( 0.09, "X");

  MT2StyleRatio->SetTitleSize  ( 0.150,"Y");
  MT2StyleRatio->SetTitleOffset( 0.51,"Y"); // used to be  0.44
  MT2StyleRatio->SetLabelSize  ( 0.150,"Y");
  MT2StyleRatio->SetNdivisions ( 505  ,"Y");
  MT2StyleRatio->SetLabelFont  ( 42   ,"Y");
  MT2StyleRatio->SetTitleFont  ( 42   ,"Y");
  MT2StyleRatio->SetTickLength (0.04,  "Y");

  MT2StyleRatio->SetTitleSize  ( 0.150,"Z");
  MT2StyleRatio->SetTitleOffset( 1.100,"Z");
  MT2StyleRatio->SetLabelOffset( 0.015,"Z");
  MT2StyleRatio->SetLabelSize  ( 0.150,"Z");
  MT2StyleRatio->SetLabelFont  ( 42   ,"Z");
  MT2StyleRatio->SetTitleFont  ( 42   ,"Z");
  MT2StyleRatio->SetNdivisions ( 707  ,"Z");


  MT2StyleRatio->SetTitleBorderSize  (0);
  MT2StyleRatio->SetTitleFillColor  (0);  
  MT2StyleRatio->SetTitleFont  (42);
  MT2StyleRatio->SetTitleColor  (1);
  MT2StyleRatio->SetLineWidth  (1);

  // Options
  MT2StyleRatio->SetOptFit        (0);
  MT2StyleRatio->SetOptStat       (0);
  MT2StyleRatio->SetOptTitle      (0);
  MT2StyleRatio->SetStatBorderSize(0);
  MT2StyleRatio->SetStatColor     (0);
  MT2StyleRatio->SetStatW         (0.5);

// Default palette gradient - from Marco-Andrea / Marco Rossini
  Double_t stops[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red  [5] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[5] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue [5] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(5, stops, red, green, blue, 100); 
  MT2StyleRatio->SetNumberContours(100);

}
