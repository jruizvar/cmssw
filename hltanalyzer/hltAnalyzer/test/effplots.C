// Usage
// root -b -q 'effplots.C("eeMass")'
//
void effplots(char* plot)
{
gROOT->SetStyle("Plain");
TFile *f1 = new TFile("histoGenElectrons_M1000PU20bx25.root");
TFile *f2 = new TFile("histoGenElectrons_M1000PU40bx25.root");
TFile *f3 = new TFile("histoGenElectrons_M1000PU40bx50.root");
TEfficiency *e1, *e2, *e3;
f1->GetObject(Form("demo/eff_%s_HLT_DoubleEle33_CaloIdL_v15;1",plot),e1);
f2->GetObject(Form("demo/eff_%s_HLT_DoubleEle33_CaloIdL_v15;1",plot),e2);
f3->GetObject(Form("demo/eff_%s_HLT_DoubleEle33_CaloIdL_v15;1",plot),e3);
e1->SetLineColor(4); e1->SetLineWidth(1); e1->SetMarkerStyle(24); e1->SetMarkerSize(0.8); e1->SetMarkerColor(4); 
e2->SetLineColor(2); e2->SetLineWidth(2); e2->SetMarkerStyle(25); e2->SetMarkerSize(0.8); e2->SetMarkerColor(2);    
e3->SetLineColor(1); e3->SetLineWidth(2); e3->SetMarkerStyle(26); e3->SetMarkerSize(0.8); e3->SetMarkerColor(1); 
c1 = new TCanvas("c1");
c1->cd();
e1->Draw();
gPad->SetGridx();
gPad->Update();
e1->GetPaintedGraph()->SetMaximum(1.);
e1->GetPaintedGraph()->SetMinimum(0.);
e2->Draw("Xsame");
e3->Draw("Xsame");
leg = new TLegend(0.83,0.77,0.99,0.99);
leg->SetTextSize(0.035);
leg->SetFillColor(0);
leg->SetLineColor(1);
leg->SetHeader("RSGrav M1000");
leg->AddEntry(e1,"PU20bx25","lp");
leg->AddEntry(e2,"PU40bx25","p");
leg->AddEntry(e3,"PU40bx50","p");
leg->Draw();
c1->Print(Form("eff_%s_M1000.png",plot));

f1 = new TFile("histoGenElectrons_M4500PU20bx25.root");
f2 = new TFile("histoGenElectrons_M4500PU40bx25.root");
f3 = new TFile("histoGenElectrons_M4500PU40bx50.root");
f1->GetObject(Form("demo/eff_%s_HLT_DoubleEle33_CaloIdL_v15;1",plot),e1);
f2->GetObject(Form("demo/eff_%s_HLT_DoubleEle33_CaloIdL_v15;1",plot),e2);
f3->GetObject(Form("demo/eff_%s_HLT_DoubleEle33_CaloIdL_v15;1",plot),e3);
e1->SetLineColor(4); e1->SetLineWidth(1); e1->SetMarkerStyle(24); e1->SetMarkerSize(0.8); e1->SetMarkerColor(4); 
e2->SetLineColor(2); e2->SetLineWidth(2); e2->SetMarkerStyle(25); e2->SetMarkerSize(0.8); e2->SetMarkerColor(2);    
e3->SetLineColor(1); e3->SetLineWidth(2); e3->SetMarkerStyle(26); e3->SetMarkerSize(0.8); e3->SetMarkerColor(1); 
c2 = new TCanvas("c2","");
c2->cd();
e1->Draw();
gPad->SetGridx();
gPad->Update();
e1->GetPaintedGraph()->SetMaximum(1.);
e1->GetPaintedGraph()->SetMinimum(0.);
e2->Draw("Xsame");
e3->Draw("Xsame");
leg->SetHeader("RSGrav M4500");
leg->Draw();
c2->Print(Form("eff_%s_M4500.png",plot));
}
