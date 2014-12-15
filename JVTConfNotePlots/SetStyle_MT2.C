void SetStyle_MT2(){

  TStyle *RootStyle = new TStyle("MT2-Style","MT2-Style");
  
#ifdef __CINT__
  TStyle *GloStyle;
  GloStyle = gStyle;                          // save the global style reference
  gStyle = RootStyle;
#endif


  // Canvas

  RootStyle->SetCanvasColor     (0);
  RootStyle->SetCanvasBorderMode(0);
  RootStyle->SetCanvasBorderMode(0);
  RootStyle->SetCanvasDefH      (600);
  RootStyle->SetCanvasDefW      (600);

  // Pads

  RootStyle->SetPadColor       (0);
  RootStyle->SetPadBorderMode  (0);
  RootStyle->SetPadBorderSize  (1);
  RootStyle->SetPadBottomMargin(0.17);
  RootStyle->SetPadTopMargin   (0.07);
  RootStyle->SetPadLeftMargin  (0.18);
  RootStyle->SetPadRightMargin (0.08);
  RootStyle->SetPadTickX       (1);
  RootStyle->SetPadTickY       (1);

  // Frames
  RootStyle->SetFrameLineWidth (1);
  RootStyle->SetFrameFillColor (0);
  RootStyle->SetFrameBorderMode(0);
  RootStyle->SetFrameBorderSize(0);


  // Histograms
   RootStyle->SetHistFillColor(0);
   RootStyle->SetHistLineWidth(1);

  // Functions
  //RootStyle->SetFuncColor(1);
  //RootStyle->SetFuncStyle(0);
  //RootStyle->SetFuncWidth(2);

  //Legends 
  RootStyle->SetLegendBorderSize(0);
  RootStyle->SetFillStyle(0);
  RootStyle->SetLegendFont(42); 

  // Labels, Ticks, and Titles
  RootStyle->SetTitleSize  ( 0.050,"X");
  RootStyle->SetTitleOffset( 1.100,"X");
  RootStyle->SetLabelSize  ( 0.050,"X");
  RootStyle->SetLabelFont  ( 42   ,"X");
  RootStyle->SetTitleFont  ( 42   ,"X");
  RootStyle->SetNdivisions ( 505  ,"X");

  RootStyle->SetTitleSize  ( 0.050,"Y");
  RootStyle->SetTitleOffset( 1.500,"Y");
  RootStyle->SetLabelSize  ( 0.050,"Y");
  RootStyle->SetNdivisions ( 505  ,"Y");
  RootStyle->SetLabelFont  ( 42   ,"Y");
  RootStyle->SetTitleFont  ( 42   ,"Y");

  RootStyle->SetTitleSize  ( 0.050,"Z");
  RootStyle->SetTitleOffset( 1.100,"Z");
  RootStyle->SetLabelOffset( 0.015,"Z");
  RootStyle->SetLabelSize  ( 0.050,"Z");
  RootStyle->SetLabelFont  ( 42   ,"Z");
  RootStyle->SetTitleFont  ( 42   ,"Z");
  RootStyle->SetNdivisions ( 707  ,"Z");


  RootStyle->SetTitleBorderSize  (0);
  RootStyle->SetTitleFillColor  (0);  
  RootStyle->SetTitleFont  (42);
  RootStyle->SetTitleColor  (1);

  RootStyle->SetLineWidth  (1);

  // Options
  RootStyle->SetOptFit        (0);
  RootStyle->SetOptStat       (0);
  RootStyle->SetOptTitle      (0);
  RootStyle->SetStatBorderSize(0);
  RootStyle->SetStatColor     (0);
  RootStyle->SetStatW         (0.5);

 ////////////////// Default palette gradient - from Marco-Andrea / Marco Rossini
  Double_t stops[5] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red  [5] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[5] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue [5] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(5, stops, red, green, blue, 100); 
  RootStyle->SetNumberContours(100);

}
