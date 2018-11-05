#include<iostream>
#include<vector>
#include<algorithm>
#include<TH1F.h>
#include<TCanvas.h>
#include<map>

void merge(){

  gROOT->SetBatch(kTRUE);

  std::cout<<"so fucked up"<<std::endl;  

  std::vector<int> key_bottom { 11,12,14,17,18,19,22,23,24 }; 
  std::vector<int> key_ft     { 26,27,28,29,30,31,52,56,57,58,59,60,61 }; 
  std::vector<int> key_pipe   { 15,16,20,21,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,53,54,55 }; 
  std::vector<int> key_top    { 105,106,107,108,109,111,112,113,114,115,116,117,118,119,120,121,195,123,124,125,126,127,128,129 }; //Top 122 is using FEB 195

  std::vector<std::vector<int>> key_planes {key_bottom, key_ft, key_pipe, key_top};
  std::vector<std::string> path {"bottom", "ft", "pipe","top"};
  
  
  TFile *_file_data = TFile::Open("strip_tree_data.root");
  TFile *_file_mc   = TFile::Open("strip_tree_mc.root");
  
  int n_hist_generated = 0 ;

  for (int plane=0; plane<4; ++plane){
    int key_feb=0;
    //Bottom
    for(auto feb_id : key_planes[plane]){
      key_feb=feb_id*100;
      for (int strip=0; strip<16; ++strip){
	int strip_id=key_feb+strip;
	n_hist_generated++;
	std::string s = std::to_string(strip_id);
	char const *id_char = s.c_str();

	std::cout<<strip_id<<std::endl;
	
	for (int ax=0; ax<3; ++ax){
	  char ax_label[] = "_x";
	  if(ax==1) ax_label[1] = 'y';
	  if(ax==2) ax_label[1] = 'z';
	  char * title = new char[std::strlen(id_char)+std::strlen(ax_label)+1];
	  std::strcpy(title,id_char);
	  std::strcat(title,ax_label);
	  


	  TH1F * hist_data= (TH1F*)   _file_data->Get(Form("crt/%s/%s",path[plane].c_str(),title));
	  TH1F * hist_mc  = (TH1F*)   _file_mc  ->Get(Form("crt/%s/%s",path[plane].c_str(),title));
	  
	  TCanvas *c_out = new TCanvas(title,"strip",500,500);
	  c_out->cd();
	  hist_data->Draw();
	  hist_mc->SetLineColor(kRed);
	  hist_mc->Draw("same");

	  auto legend = new TLegend(0.79,0.64,0.89,0.74);
	  //legend->SetHeader(title,"C");
	  legend->AddEntry(hist_data,"Data","l");
	  legend->AddEntry(hist_mc  ,"MC"  ,"l");
	  legend->Draw();
	  
	  c_out->SaveAs(Form("%s/%s.pdf", path[plane].c_str(),title));
	}//per ub plane
      }//per strip
    }//per feb
  }//per CRT plane
}
