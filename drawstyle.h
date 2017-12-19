void setPad(Double_t left, Double_t right, Double_t top, Double_t bottom, int color=10)
{
    gPad->SetFillColor(color);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFrameFillColor(10);
    gPad->SetFrameBorderMode(0);
    gPad->SetFrameBorderSize(0);
    gPad->SetLeftMargin(left);
    gPad->SetRightMargin(right);
    gPad->SetTopMargin(top);
    gPad->SetBottomMargin(bottom);
}

void setGraphStyle(TGraphErrors* gr, int color, int marker, int mSize, int lSize)
{
    gr->SetLineColor(color);
    gr->SetMarkerColor(color);
    gr->SetMarkerStyle(marker);
    gr->SetMarkerSize(mSize);
    gr->SetLineWidth(lSize);
}

//__________________________________________________
TLine* drawLine(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineStyle,Int_t lineColor)
{
    TLine *l1 = new TLine(xlow,ylow,xup,yup);
    l1->SetLineWidth(lineWidth);
    l1->SetLineColor(lineColor);
    l1->SetLineStyle(lineStyle);
    l1->Draw("same");
    return l1;
}

void myStyle() {
    
        cout << "Welcome to Shuai's style Setting" << endl;
    
        TStyle *myStyle= new TStyle("myStyle","my plots style");
    
        myStyle->SetPalette(1,0);
    
        // use plain black on white colors
        myStyle->SetCanvasColor(10);
        myStyle->SetCanvasBorderMode(0);
        myStyle->SetCanvasBorderSize(2);
        myStyle->SetPadColor(10);
        myStyle->SetPadBorderMode(0);
        myStyle->SetPadBorderSize(0);
        myStyle->SetPadBottomMargin(0.12);
        myStyle->SetPadLeftMargin(0.14);
        myStyle->SetPadRightMargin(0.08);
        myStyle->SetPadTopMargin(0.08);
        myStyle->SetLineWidth(2);//change tick width
        myStyle->SetPadTickX(1);
        myStyle->SetPadTickY(1);
        myStyle->SetTickLength(0.02,"X");
        myStyle->SetTickLength(0.02,"Y");
        myStyle->SetPadGridX(0);
        myStyle->SetPadGridY(0);
        myStyle->SetGridColor(18);
        myStyle->SetFrameFillStyle(4000);
        myStyle->SetFrameLineWidth(2);
        myStyle->SetFrameBorderSize(2);
        myStyle->SetFrameBorderMode(0);
        myStyle->SetFrameFillColor(10);
        //gStyle->SetFrameLineStyle(1);
        myStyle->SetLegendBorderSize(0);
    
        // set the paper & margin sizes
        myStyle->SetPaperSize(20,26);
    
        int font = 22;
        // use large Times-Roman fonts
        myStyle->SetTextFont(font);
        myStyle->SetTextSize(0.08);
        myStyle->SetLabelFont(font,"xyz");
        myStyle->SetTitleFont(font,"xyz");
        myStyle->SetLegendFont(font);
        myStyle->SetStatFont(font);
    
        myStyle->SetLabelSize(0.04,"xyz");//D=0.04
        myStyle->SetTitleSize(0.05,"xyz");//D=0.02
        myStyle->SetTitleSize(0.05,"");//main title
    
        myStyle->SetLabelOffset(0.015,"xyz");//D=0.005
        myStyle->SetTitleOffset(1.0,"x");
        myStyle->SetTitleOffset(1.2,"y");
        myStyle->SetTitleOffset(1.2,"z");
        TGaxis::SetMaxDigits(3);
    
        // use bold lines and markers
        myStyle->SetMarkerStyle(20);
        myStyle->SetHistLineWidth(1.2);
        myStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
    
        // get rid of X error bars and y error bar caps
        //myStyle->SetErrorX(0.001);
    
        // do not display any of the standard histogram decorations
        myStyle->SetTitleX(0.5);
        myStyle->SetTitleAlign(23);
        //myStyle->SetTitleColor(0);
        myStyle->SetTitleStyle(0);
        myStyle->SetTitleBorderSize(0);
        myStyle->SetOptTitle(0);
        myStyle->SetOptStat(0);
        myStyle->SetOptFit(0);
    
        myStyle->SetStatColor(10);
    
        gROOT->SetStyle("myStyle");
        gROOT->ForceStyle();
    
    }
    
    void globleStyle() {
    
        //============================================================
        // 
        //          Make graphs pretty
        // 
        //============================================================
        cout<<"gStyle mode requested!!!"<<endl;
    
        int font = 22;
    
        gStyle->SetOptTitle(1);
        gStyle->SetOptDate(0);
        gStyle->SetOptStat(0);
        gStyle->SetStatColor(10);
        //gStyle->SetOptFit(0);
        gStyle->SetStatH(0.17);
        gStyle->SetStatW(0.17);
        gStyle->SetPalette(1,0);
        gStyle->SetTextFont(font);
        gStyle->SetTextSize(0.055);
        //gStyle->SetErrorX(1);
        gStyle->SetEndErrorSize(4);
        gStyle->SetDrawBorder(0);
    
        gStyle->SetCanvasDefH(600);
        gStyle->SetCanvasDefW(800);
        gStyle->SetCanvasColor(10);
        gStyle->SetCanvasBorderMode(0);
        gStyle->SetCanvasBorderSize(2);
        gStyle->SetPadColor(10);
        gStyle->SetPadBorderMode(0);
        gStyle->SetPadBorderSize(0);
        gStyle->SetPadBottomMargin(0.12);
        gStyle->SetPadLeftMargin(0.12);
        gStyle->SetPadRightMargin(0.10);
        gStyle->SetPadTopMargin(0.08);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetTickLength(0.02,"X");
        gStyle->SetTickLength(0.02,"Y");
        gStyle->SetPadGridX(0);
        gStyle->SetPadGridY(0);
        gStyle->SetGridColor(18);
        gStyle->SetLineWidth(2);
        gStyle->SetFrameFillStyle(4000);
        gStyle->SetFrameLineWidth(2);
        gStyle->SetFrameBorderSize(2);
        gStyle->SetFrameBorderMode(0);
        gStyle->SetFrameFillColor(10);
        //gStyle->SetFrameLineStyle(1);
    
        gStyle->SetNdivisions(510,"X");
        gStyle->SetNdivisions(510,"Y");
        gStyle->SetLabelSize(0.04,"X");
        gStyle->SetLabelSize(0.04,"Y");
        gStyle->SetLabelFont(font,"X");
        gStyle->SetLabelFont(font,"Y");
        gStyle->SetLabelOffset(0.01,"X");
        gStyle->SetLabelOffset(0.01,"Y");
        gStyle->SetTitleOffset(1.0,"X");
        gStyle->SetTitleOffset(1.2,"Y");
        gStyle->SetTitleOffset(1.2,"Z");
        gStyle->SetTitleSize(0.05,"X");
        gStyle->SetTitleSize(0.05,"Y");
        gStyle->SetTitleSize(0.05,"Z");
        gStyle->SetTitleFont(font,"X");
        gStyle->SetTitleFont(font,"Y");
        gStyle->SetTitleFont(font,"Z");
    
        // COPY FROM MYSTYLE()
        // gStyle->SetTitleColor(3);
        // gStyle->SetTitleX(0.07);
        // gStyle->SetTitleAlign(23);
        // gStyle->SetTitleStyle(0);
        // gStyle->SetTitleBorderSize(0);
    }
    
    // set color display, raibow, grayscale...
    void set_color_env(){
    
        const Int_t NRGBs = 5;
        const Int_t NCont = 255;
    
        Int_t fcol;
    
        Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
        Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
        Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
        Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
        fcol = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
        gStyle->SetNumberContours(NCont);
        //SetPalette has been called in the above function, color style will
        //be set according to your own rgb definition if called
    
        //grayscale
        /*
           double dcol = 1/double(NRGBs);
           double grey = 1;
    
           for(int j = 0; j < NRGBs; j++){  
        // ...... Define color with RGB equal to : gray, gray, gray .......
        stops[j]=double(j)/double(NRGBs-1);
        red[j]=grey;
        blue[j]=grey;
        green[j]=grey;
    
        TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
        */	
    
        TString Red1 =  "#FFCC00";
        TString Red2 =  "#FF9900";
        TString Red3 =  "#FF6600";
        TString Red4 =  "#FF3300";
        TString Red5 =  "#FF0000";
    
        TString Blue1 =  "#3300FF";
        TString Blue2 =  "#0000FF";
        TString Blue3 =  "#0033FF";
        TString Blue4 =  "#0066FF";
        TString Blue5 =  "#0099FF";
        TString Blue6 =  "#00CCFF";
        TString Blue7 =  "#00FFFF";
        TString Blue8 =  "#00FFCC";
    
        TString Green1 =  "#006633";
        TString Green2 =  "#006600";
        TString Green3 =  "#009933";
        TString Green4 =  "#009900";
        TString Green5 =  "#339900";
        TString Green6 =  "#00CC33";
        TString Green7 =  "#00CC00";
        TString Green8 =  "#33CC00";
        TString Green9 =  "#00FF00";
    
        TString Oran1 = "#FF3300";
        TString Oran2 = "#FF6600";
        TString Oran3 = "#FF6633";
        TString Oran4 = "#FF9900";
        TString Oran5 = "#FF9933";
        TString Oran6 = "#FF9966";  
    
        TString SkyBlue = "#00CCFF";
        TString SeaBlue = "#0099FF";
        TString SadBlue = "#009999";
        TString LakeBlue = "#0099CC";
        TString DarkBlue = "#000099";
    
        TString Purp1 = "#CC33CC";
        TString Purp2 = "#9900FF";
        TString Purp3 = "#CC00FF";
        TString Purp4 = "#FF00FF";
        TString Purp5 = "#FF33FF";
        TString Purp6 = "#FF33CC";   
        TString Purp7 = "#FF66FF";
        //how to use the color defined above --->  SetLineColor(TColor::GetColor(SkyBlue)) or SetLineColor(TColor::GetColor(SkyBlue.Data()))
    
    }