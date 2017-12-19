#include "../anaconst.h"

void mixedbkgd(int lcent = 7, int hcent = 8){

	double eff_pt[nPtBins]={0};
	double sigma_pt[nPtBins]={0};

	//float const PtEdge[nPtBins+1] = {1.5, 2.,3.,4.,5.,6,8,12};
	double PtCenter[nPtBins];
	
	TString rootname = "../D0withrapiditybin.root";
	TFile *f = TFile::Open(rootname);

	TH3D* imp = NULL;
	TH3D* implk = NULL;
	
	imp= (TH3D*)f->Get("mh3InvariantMassVsPtVsCentTof")->Clone("mh3InvariantMassVsPtVsCentTof");
	implk= (TH3D*)f->Get("mh3InvariantMassVsPtVsCentTofLike")->Clone("mh3InvariantMassVsPtVsCentToflike");
	//TH1F* ref = (TH1F*)f->Get("mh1CentWg")->Clone("mh1CentWg");
	//ref->SetDirectory(0);	
	cout<<"ok2"<<endl;
	imp->SetDirectory(0);
	implk->SetDirectory(0);
	f->Close();
	
	TString bkname = "Mixed_Evt.2017Aug16.root";
	//TString bkname = "D0_Data_Mix.root";
	//TString bkname = "D0_Data_Mix_MergeUsed.root";
	f=TFile::Open(bkname);
	THnF* hmix = (THnF*)f->Get("hD0CentPtEtaMDphiDaugMixed_standard;1")->Clone("hD0CentPtEtaMDphiDaugMixed_standard");
	cout<<"ok3"<<endl;
	THnF* hmixlk = (THnF*)f->Get("hD0CentPtEtaMDphiDaugLikeSignMixed_standard")->Clone("hD0CentPtEtaMDphiDaugLikeSignMixed_standard");
	cout<<"ok4"<<endl;
	THnF* hmass = (THnF*)f->Get("hD0CentPtEtaMDphiDaug_standard;1")->Clone("hD0CentPtEtaMDphiDaug_standard");
	cout<<"ok5"<<endl;
	THnF* hmasslk = (THnF*)f->Get("hD0CentPtEtaMDphiDaugLikeSign_standard;1")->Clone("hD0CentPtEtaMDphiDaugLikeSign_standard");
	//hmix->SetDirectory(0);
	//hmixlk->SetDirectory(0);
	//hmass->SetDirectory(0);
	//hmasslk->SetDirectory(0);
	TH1F* ref = (TH1F*)f->Get("hCentralityWeighted")->Clone("mh1CentWg");
	ref->SetDirectory(0);
	f->Close();
	
	fstream ftot;
	ftot.open(Form("data/D0tot_%d_%d.txt",centedge[hcent+1],centedge[lcent]), ios::out | ios::trunc);
	long events = ref->Integral(ref->FindBin(lcent),ref->FindBin(hcent));
	ftot<<events<<endl;
	ftot.close();
	
	int cbinh = hmix->GetAxis(0)->FindBin(hcent+0.5);
	int cbinl = hmix->GetAxis(0)->FindBin(lcent+0.5);
	cout<<"cent"<<cbinh<<"\t"<<cbinl<<endl;
	hmix->GetAxis(0)->SetRange(cbinl,cbinh);  //centrality and it add 0.5 when fill it 
	hmixlk->GetAxis(0)->SetRange(cbinl,cbinh);  //centrality and it add 0.5 when fill it 
	hmass->GetAxis(0)->SetRange(cbinl,cbinh);  //centrality and it add 0.5 when fill it 
	hmasslk->GetAxis(0)->SetRange(cbinl,cbinh);  //centrality and it add 0.5 when fill it 
	
	int mbinh = hmix->GetAxis(3)->FindBin(2.1-1e-6);
	int mbinl = hmix->GetAxis(3)->FindBin(1.6+1e-6);
	hmix->GetAxis(3)->SetRange(mbinl,mbinh);
	hmixlk->GetAxis(3)->SetRange(mbinl,mbinh);
	hmass->GetAxis(3)->SetRange(mbinl,mbinh);
	hmasslk->GetAxis(3)->SetRange(mbinl,mbinh);
	
	double ptmin = 0.3;
	int dau1PtBin = hmix->GetAxis(2)->FindBin(ptmin+1e-6);
	int dau2PtBin = hmix->GetAxis(4)->FindBin(ptmin+1e-6);
	int dau1NPtBin = 6;  //11
	int dau2NPtBin = 6;  //11
	hmix->GetAxis(2)->SetRange( dau1PtBin, dau1NPtBin);  //
	hmix->GetAxis(4)->SetRange( dau2PtBin, dau2NPtBin);  //
	hmixlk->GetAxis(2)->SetRange( dau1PtBin, dau1NPtBin);  //
	hmixlk->GetAxis(4)->SetRange( dau2PtBin, dau2NPtBin);  //
	hmass->GetAxis(2)->SetRange( dau1PtBin, dau1NPtBin);  //
	hmass->GetAxis(4)->SetRange( dau2PtBin, dau2NPtBin);  //
	hmasslk->GetAxis(2)->SetRange( dau1PtBin, dau1NPtBin);  //
	hmasslk->GetAxis(4)->SetRange( dau2PtBin, dau2NPtBin);  //
	cout<<"ok"<<endl;
	for (int i=0;i<nPtBins;i++){
		int ptbinh = hmix->GetAxis(1)->FindBin(PtEdge[i+1]-1e-6);
		int ptbinl = hmix->GetAxis(1)->FindBin(PtEdge[i]+1e-6);

		hmass->GetAxis(1)->SetRange( ptbinl, ptbinh);  //
		hmasslk->GetAxis(1)->SetRange( ptbinl, ptbinh);  //
		hmix->GetAxis(1)->SetRange( ptbinl, ptbinh);  //
		hmixlk->GetAxis(1)->SetRange( ptbinl, ptbinh);  //
		TH1D* h = imp->ProjectionZ("imp_pz",imp->GetXaxis()->FindBin(PtEdge[i]+1e-6),  imp->GetXaxis()->FindBin(PtEdge[i+1]-1e-6), imp->GetYaxis()->FindBin(lcent),  imp->GetYaxis()->FindBin(hcent),"e");
		TH1D* hlk = implk->ProjectionZ("implk_pz",imp->GetXaxis()->FindBin(PtEdge[i]+1e-6),  imp->GetXaxis()->FindBin(PtEdge[i+1]-1e-6), imp->GetYaxis()->FindBin(lcent),  imp->GetYaxis()->FindBin(hcent),"e");
		//TH1D* h = hmass->Projection(3,"E");
		//TH1D* hlk = hmasslk->Projection(3,"E");
		
		TH1D* hmixbk = hmix->Projection(3,"E");
		TH1D* hmixbklk = hmixlk->Projection(3,"E");
		hmixbk->SetLineColor(kBlue);
		
		double scale=1;
		std::pair<float, float> sideband1(1.6, 1.8);
		std::pair<float, float> sideband2(1.93, 2.1);
		std::pair<float, float> sideband(1.6, 2.1);
//		double sidemix = hmixbk->Integral(hmixbk->FindBin(sideband1.first), hmixbk->FindBin(sideband1.second)) +
//								  hmixbk->Integral(hmixbk->FindBin(sideband2.first), hmixbk->FindBin(sideband2.second));
//		double sidesignal = h->Integral(h->FindBin(sideband1.first),h->FindBin(sideband1.second)) +
//									 h->Integral(h->FindBin(sideband2.first),h->FindBin(sideband2.second));
		double sidemix = hmixbk->Integral(hmixbk->FindBin(sideband.first+1.e-6), hmixbk->FindBin(sideband.second-1.e-6));
		double sidemixlk = hmixbklk->Integral(hmixbklk->FindBin(sideband.first+1.e-6), hmixbk->FindBin(sideband.second-1.e-6));
		double sidesignal = hlk->Integral(hlk->FindBin(sideband.first+1.e-6),hlk->FindBin(sideband.second-1e-6));
		//scale = sidemix? sidesignal/sidemix;
		scale = sidemixlk? sidesignal/sidemixlk;
		cout<<scale<<endl;
		hmixbk->Scale(scale);
		
		TH1F* hdelta = (TH1F*)h->Clone("delta");
		hdelta->Sumw2();
		hdelta->Add(h, hmixbk, 1, -1);
		hlk->Draw();
		h->Draw("same");
		h->SetLineColor(kRed);
		h->SetMarkerColor(kRed);
		hmixbk->Draw("same");
		hdelta->Draw("same");
		fit(hdelta,sigma_pt[i],eff_pt[i], lcent, hcent);
		PtCenter[i]=(PtEdge[i]+PtEdge[i+1])/2.0;                                                      
	}
	
	TCanvas* c = new TCanvas("cc","cc");
	TFile * fout = new TFile("D0masscuteff_mix.root","recreate");
	TH1D* hmasseff= new TH1D("hmasseff","masseff;p_{T}/GeV/c;Eff.", nPtBins, PtEdge); 
	for (int i=0;i<nPtBins;i++){
		hmasseff->SetBinContent(i+1,eff_pt[i]);
		hmasseff->SetBinError(i+1,sigma_pt[i]);
	}
	
	hmasseff->Draw();
	hmasseff->Write();
	
	fout->Close();
	gPad->Close();
}

Double_t  myGaus(Double_t *x,Double_t *par){ 
	Double_t fitval=0;
	Double_t arg = 0;
	if (par[2]!=0) arg = (x[0] - par[1])/par[2];
	Double_t bin_width = 0.01;
	//Double_t constant = bin_width*par[0]/sqrt(2*TMath::Pi())/par[2];
	//fitval = constant*TMath::Exp(-0.5*arg*arg)+par[3]+par[4]*TMath::Exp(par[5]*x[0]);
	//fitval = bin_width*par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE)+par[3]+par[4]*TMath::Exp(par[5]*x[0]);
	fitval = bin_width*par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE)+par[3]+par[4]*x[0];
	return fitval;
}

std::pair<float, float> gMassExclusionRange(1.82, 1.9);
Double_t funResidualBg(Double_t *x, Double_t *par)
{
    if (x[0] > gMassExclusionRange.first && x[0] < gMassExclusionRange.second)
    {
        TF1::RejectPoint();
        return 0;
    }
    
    return par[0] + par[1] * x[0];
}

void SetHist(TH3* h){
	h->GetXaxis()->SetTitle("p_{T}(GeV/c)");
	h->GetYaxis()->SetTitle("Cent");
	h->GetZaxis()->SetTitle("m_{K#pi}(GeV/c^{2})");
	h->Sumw2();
}

void fit(TH1* h,double &sigma,double &eff, int lcent = 7, int hcent = 8){
	double fitrang_low = 1.72;
	double fitrange_up = 2.1;
	h->Draw();
	TF1* gaus = new TF1("gaus","gaus",1.83,1.9);
	h->Fit(gaus,"R");
	const int npar = 5;
	double par[npar]={0};
	gaus->GetParameters(par);
	// TF1*  bk= new TF1("bk","[0]+[1]*TMath::Exp([2]*x)",1.75,2.1);
    TF1*  bk= new TF1("bk", funResidualBg,1.6,2.1, 2);
//	TF1*  bk= new TF1("bk", "pol1(0)",1.91,2.1);
	bk->SetParameters(0,0);
	h->Fit(bk,"IONR","",fitrang_low, fitrange_up);
	h->Fit(bk,"IONR","",fitrang_low, fitrange_up);
	bk->GetParameters(par+3);
	TF1* mygaus=new TF1("mygaus",myGaus,1.6,2.1,npar);
	//TF1* mygaus=new TF1("mygaus","0.01*[0]*gausn(1)+pol1(4)",1.6,2.1);
	mygaus->SetParameters(par);
	h->Fit(mygaus,"IOR", "", fitrang_low, fitrange_up);
	//sigma=mygaus->GetParameter(2);
	mygaus->GetParameters(par);
	bk->SetParameters(par+3);
	bk->Draw("same");
	TF1* normgaus = new TF1("normgaus", "0.01*[0]*TMath::Gaus(x,[1],[2],1)",  1.65, 2.1 );
	normgaus->SetParameters(par);
	normgaus->Draw("same");
	normgaus->SetLineColor(kGreen);
	//double tot=normgaus->Integral(1.65,2.1);
	//double tot=h->Integral(h->FindBin(1.75),h->FindBin(2.1))-bk->Integral(1.75,2.1)/h->GetBinWidth(1);
	//double part=h->Integral(h->FindBin(1.83),h->FindBin(1.9))-bk->Integral(1.83,1.9)/h->GetBinWidth(1);
	double tot=mygaus->GetParameter(0);
	ofstream out;
	out.open(Form("data/D0tot_%d_%d.txt",centedge[hcent+1],centedge[lcent]),ios::app);
	out<<tot<<"\t"<<mygaus->GetParError(0)<<endl;
	double part=normgaus->Integral(1.83,1.90);
	eff=(part+1)/(tot+2);
	sigma = sqrt(eff*(tot-part+1)/(tot+2)/(tot+3));
	gPad->Update();
	gSystem->Sleep(500);
	delete gaus;
	delete bk;
	delete mygaus;
}

