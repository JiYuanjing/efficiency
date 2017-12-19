#include "anaconst.h"

void IntRatio(){
 fstream fratio;
 fratio.open("data/intratio.txt", ios::trunc | ios::out );
 for (int i=0;i<nCentBins;i++)
 {
  int hcent = hCentEdge[i];
  int lcent = lCentEdge[i];
  cout<<centedge[hcent+1]<<"\t"<<centedge[lcent]<<endl;
  TFile* fD0 = TFile::Open(Form("D0/data/hD0_myspectra_%d_%d.root", centedge[hcent+1],centedge[lcent]));
  TH1D* hD0 = (TH1D*)fD0->Get("D0effyield")->Clone(Form("D0effyield_%d_%d",centedge[hcent+1],centedge[lcent]));
  hD0->SetDirectory(0);
  fD0->Close();
  TFile* fDstar = TFile::Open(Form("data/correct_%d_%d.root",centedge[hcent+1],centedge[lcent]));
  TH1D* hDstar = (TH1D*)fDstar->Get(Form("heffyield_%d_%d",centedge[hcent+1],centedge[lcent]))->Clone(Form("heffyield_%d_%d",centedge[hcent+1],centedge[lcent]));
  hDstar->SetDirectory(0);
  fDstar->Close();
  double hpt = 10,lpt=2;
  double NDstar, NDstarerr, ND0, ND0err;
  CountsandErr(hpt, lpt, ND0, ND0err, hD0);
  CountsandErr(hpt, lpt, NDstar, NDstarerr, hDstar);
  cout<<ND0<<"\t"<<ND0err<<"\t"<<endl;
  double ratio = NDstar/ND0/DstarBranch;
  double error = sqrt(1/ND0err/ND0err+NDstarerr*NDstarerr/ND0err/ND0err/ND0err/ND0err)/DstarBranch;
  fratio<<centedge[hcent+1]<<"\t"<<centedge[lcent]<<"\t"<<lpt<<"\t"<<hpt<<"\t"<<ratio<<"\t"<<error<<endl;
 }
 fratio<<"minCent"<<"\t"<<"maxCent"<<"\t"<<"lowPt"<<"\t"<<"highPt"<<"\t"<<"ratio"<<endl;
 fratio.close();
}

void CountsandErr(double hpt, double lpt, double& N, double& Nerr, TH1* h ){
  N = h->Integral(h->FindBin(lpt+1e-6),h->FindBin(hpt-1e-6));
  for (int i=h->FindBin(lpt+1e-6);i<=h->FindBin(hpt-1e-6);i++){
    Nerr+=h->GetBinError(i);
  }
  Nerr = sqrt(Nerr);
}
