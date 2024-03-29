/////////
// the Style Section
/////////


/////////
// Colors, markers, and line styles that are easy to differentiate in 
// slides, and color / black and white print
/////////

//  int color[] = { 1, 2, 4, 904, 419, 9,11,12,13,14,15,16,17,18,19,20,21};
//  int style[] = {21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38};
//  int lstyle[] ={ 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9};


void setTDRStyle() {
  //  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  cout << "Setting TDRStyle" << endl;

// For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(500); //Height of canvas for form=1 //600
  gStyle->SetCanvasDefW(700); //Width of canvas  for form=1 //600
  gStyle->SetCanvasDefX(0);   //POsition on screen
  gStyle->SetCanvasDefY(0);

// For the Pad:
  gStyle->SetPadBorderMode(0);
  // gStyle->SetPadBorderSize(Width_t size = 1);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

// For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);

// For the histo:
  // gStyle->SetHistFillColor(1);
  // gStyle->SetHistFillStyle(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);
  // gStyle->SetLegoInnerR(Float_t rad = 0.5);
  // gStyle->SetNumberContours(Int_t number = 20);

  gStyle->SetEndErrorSize(2);
  //gStyle->SetErrorMarker(20);
  gStyle->SetErrorX(0.);
  
  gStyle->SetMarkerStyle(20);

//For the fit/function:
  gStyle->SetOptFit(1);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);

//For the date:
  gStyle->SetOptDate(0);
  // gStyle->SetDateX(Float_t x = 0.01);
  // gStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.025);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.1);
  gStyle->SetStatW(0.15);
  // gStyle->SetStatStyle(Style_t style = 1001);
  // gStyle->SetStatX(Float_t x = 0);
  // gStyle->SetStatY(Float_t y = 0);

// Margins:
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.15); //aaa 0.13
  gStyle->SetPadRightMargin(0.12); //aaa 0.05

// For the Global title:

  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
  // gStyle->SetTitleH(0); // Set the height of the title box
  // gStyle->SetTitleW(0); // Set the width of the title box
  // gStyle->SetTitleX(0); // Set the position of the title box
  // gStyle->SetTitleY(0.985); // Set the position of the title box
  // gStyle->SetTitleStyle(Style_t style = 1001);
  // gStyle->SetTitleBorderSize(2);

// For the axis titles:

  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  // gStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // gStyle->SetTitleYSize(Float_t size = 0.02);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.17); //aaa 1.05
  // gStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(true);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(true);

// Change for log plots:
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);

// Postscript options:
  // gStyle->SetPaperSize(15.,15.);
  // gStyle->SetLineScalePS(Float_t scale = 3);
  // gStyle->SetLineStyleString(Int_t i, const char* text);
  // gStyle->SetHeaderPS(const char* header);
  // gStyle->SetTitlePS(const char* pstitle);

  // gStyle->SetBarOffset(Float_t baroff = 0.5);
  // gStyle->SetBarWidth(Float_t barwidth = 0.5);
  // gStyle->SetPaintTextFormat(const char* format = "g");
  // gStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // gStyle->SetTimeOffset(Double_t toffset);
  // gStyle->SetHistMinimumZero(kTRUE);

  //  tdrStyle->cd();
  //  cout << "Finished setting TDRStyle" << endl;

}

void setRBellanStyle() {
  TStyle *theStyle = new TStyle("rbStyle", "Style for Bellan Thesis");
  theStyle->SetOptStat(0);
  theStyle->SetPadBorderMode(0);
  theStyle->SetCanvasBorderMode(0);
  theStyle->SetPadColor(0);
  theStyle->SetCanvasColor(0);
  theStyle->SetMarkerStyle(8);
  theStyle->SetMarkerSize(0.7);
  theStyle->SetPalette(1);
  
  theStyle->SetStatH(0.3);
  //   theStyle->SetTextFont(132);
  //   theStyle->SetTitleFont(132);
  theStyle->SetTitleBorderSize(1);
  //    theStyle->SetPalette(1);
  theStyle->SetOptStat(0);
  theStyle->SetFitFormat("4.4g");
  theStyle->SetStatY(0.99);
  theStyle->SetStatX(0.99);
  theStyle->SetTitleYOffset(1.6);
  theStyle->SetLabelSize(0.035, "XYZ");
  theStyle->SetPadGridX(true);
  theStyle->SetPadGridY(true);
  theStyle->SetFrameBorderMode(0);
  theStyle->SetTitleFillColor(0);
  theStyle->SetLegendBorderSize();
  
  // theStyle->SetCanvasDefH(600);
  // theStyle->SetCanvasDefW(400);
  
  //theStyle->SetOptLogy(); //aaa
  // theStyle->SetOptLogx();
  theStyle->cd();
}
