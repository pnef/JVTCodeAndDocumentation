void SetStyle_MT2_2D_newPalette(){

  TStyle *RootStyle = new TStyle("MT2_2D_newPalette","MT2_2D_newPalette");
  
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
  RootStyle->SetPadRightMargin (0.20);
  RootStyle->SetPadTickX       (0);
  RootStyle->SetPadTickY       (0);

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
  RootStyle->SetLegendFont(132); 

  // Labels, Ticks, and Titles
  RootStyle->SetTitleSize  ( 0.050,"X");
  RootStyle->SetTitleOffset( 1.100,"X");
  RootStyle->SetLabelSize  ( 0.050,"X");
  RootStyle->SetLabelFont  ( 132   ,"X");
  RootStyle->SetTitleFont  ( 132   ,"X");
  RootStyle->SetNdivisions ( 505  ,"X");

  RootStyle->SetTitleSize  ( 0.050,"Y");
  RootStyle->SetTitleOffset( 1.500,"Y");
  RootStyle->SetLabelSize  ( 0.050,"Y");
  RootStyle->SetNdivisions ( 505  ,"Y");
  RootStyle->SetLabelFont  ( 132   ,"Y");
  RootStyle->SetTitleFont  ( 132   ,"Y");

  RootStyle->SetTitleSize  ( 0.050,"Z");
  RootStyle->SetTitleOffset( 1.500,"Z");
  RootStyle->SetLabelOffset( 0.015,"Z");
  RootStyle->SetLabelSize  ( 0.050,"Z");
  RootStyle->SetLabelFont  ( 132   ,"Z");
  RootStyle->SetTitleFont  ( 132   ,"Z");
  RootStyle->SetNdivisions ( 707  ,"Z");


  RootStyle->SetTitleBorderSize  (0);
  RootStyle->SetTitleFillColor  (0);  
  RootStyle->SetTitleFont  (132);
  RootStyle->SetTitleColor  (1);

  RootStyle->SetLineWidth  (1);

  // Options
  RootStyle->SetOptFit        (0111);
  RootStyle->SetOptStat       (0);
  RootStyle->SetOptTitle      (0);
  RootStyle->SetStatBorderSize(0);
  RootStyle->SetStatColor     (0);
  RootStyle->SetStatW         (0.5);

 const int NRGBs = 4;
 const int NCont = 999;
 Double_t stops[NRGBs] = { 0.00, 0.33, 0.66, 1.00};
 Double_t red[NRGBs]   = { 1.00, 1.0,  1.0,  1.00};
 Double_t green[NRGBs] = { 1.00, 1.0,  0.16, 0.67};
 Double_t blue[NRGBs]  = { 1.00, 0.0,  0.0,  0.00};
 TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
 RootStyle->SetNumberContours(NCont);

}
