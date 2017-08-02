#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream> 
#include <cmath>
#include <exception>
#include <ctime>
#include <cmath>
#include "TNtuple.h"
#include "TROOT.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TObject.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TStyle.h"
#include "TROOT.h"
#include "THashList.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLegend.h"
#include "TProfile.h"

using namespace std;

int main(int argc, char *argv[]){

  //gStyle->SetOptStat(0);
  //gROOT->ForceStyle();
  TH1::SetDefaultSumw2();
  if(argc != 5)
    throw std::runtime_error("Bad number of arguments!");
    
  TString file1 = argv[1];
  //  TString file2 = argv[2];
  TString subDet = argv[2];
  TString slayer = argv[3];
  TString dir = argv[4];

  int layer = slayer.Atoi();

  uint16_t sdId = 0;

  if(subDet == "TIB")
    sdId = 3;
  else if(subDet == "TID")
    sdId = 4;
  else if(subDet == "TOB")
    sdId = 5;
  else if(subDet == "TEC")
    sdId = 6;
  else
    throw std::runtime_error("Wrong partition entered");
  TFile* f1 = NULL;
  TTree* t1 = NULL;
  f1 = TFile::Open(file1); 
  if(f1==NULL)
    throw std::runtime_error("File 1 address not set");
  t1 = dynamic_cast<TTree*>(f1->Get("testTree/tree"));
  if(t1==NULL)
    throw std::runtime_error("Tree 1 address not set");

  vector<float>* partition = 0;
  vector<float>* clustercharge = 0;
  vector<float>* clusterwidth = 0;
  vector<float>* clusterlayerwheel = 0;
  vector<float>* clusterstripChargeSubdetid = 0;
  vector<float>* clusterstripCharge = 0;
  vector<float>* clusterstripChargeLayerwheel = 0;
  vector<float>* clusterstripChargeStripNr = 0;
  vector<float>* clusterstripChargeTotWidth = 0;
  vector<float>* clusterstripChargeTotCharge = 0;
  vector<float>* clusterstripChargeLocalTrackPhi = 0;
  vector<float>* clusterstripChargeGlobalTrackPhi = 0;
  vector<float>* clusterstripChargeLocalTrackTheta = 0;
  vector<float>* clusterstripChargeGlobalTrackTheta = 0;
  vector<unsigned>* clusterstripChargeDetid = 0;
  vector<float>* clusterstripChargeLocalX = 0;
  vector<float>* clusterstripChargeLocalY = 0;
  vector<float>* tsostrackPt = 0;
  vector<float>* clusterstripChargetrackPt = 0;
  vector<float>* tsoslocalphi = 0;
  vector<float>* tsoslocaltheta = 0;
  vector<float>* tsosBdotY = 0;
  vector<float>* clusterstripChargelocalpitch = 0;
  vector<float>* clusterstripChargeBdotY = 0;

  vector<float> subpartition;
  vector<float> subclustercharge;
  vector<float> subclusterwidth;
  vector<float> subclusterlayerwheel;
  vector<float> subclusterstripChargeSubdetid;
  vector<float> subclusterstripCharge;
  vector<float> subclusterstripChargeLayerwheel;
  vector<float> subclusterstripChargeStripNr;
  vector<float> subclusterstripChargeTotWidth;
  vector<float> subclusterstripChargeTotCharge;
  vector<float> subclusterstripChargeLocalTrackPhi;
  vector<float> subclusterstripChargeGlobalTrackPhi;
  vector<float> subclusterstripChargeLocalTrackTheta;
  vector<float> subclusterstripChargeGlobalTrackTheta;
  vector<unsigned> subclusterstripChargeDetid;
  vector<float> subclusterstripChargeLoubclusterstripChargetrackPt;
  vector<float> subtsoslocalphi;
  vector<float> subtsoslocaltheta;
  vector<float> subtsoslocalpitch;
  vector<float> subclustersensorThickness;
  vector<float> subtsosBdotY;
  vector<float> subclusterstripChargelocalpitch;
  vector<float> subclusterstripChargesensorThickness;
  vector<float> subclusterstripChargeBdotY;



  t1->SetBranchAddress("clustersubdetid",  &partition );
  t1->SetBranchAddress("clustercharge",  &clustercharge );
  t1->SetBranchAddress("clusterwidth",  &clusterwidth );
  t1->SetBranchAddress("clusterlayerwheel",  &clusterlayerwheel );
  t1->SetBranchAddress("clusterstripChargeSubdetid",  &clusterstripChargeSubdetid );
  t1->SetBranchAddress("clusterstripCharge",  &clusterstripCharge );
  t1->SetBranchAddress("clusterstripChargeLayerwheel",  &clusterstripChargeLayerwheel );
  t1->SetBranchAddress("tsostrackPt",  &tsostrackPt );
  t1->SetBranchAddress("tsoslocalphi",  &tsoslocalphi );
  t1->SetBranchAddress("tsoslocaltheta",  &tsoslocaltheta );
  t1->SetBranchAddress("tsosBdotY",  &tsosBdotY );
  t1->SetBranchAddress("clusterstripChargeStripNr",  &clusterstripChargeStripNr );
  t1->SetBranchAddress("clusterstripChargeTotWidth",  &clusterstripChargeTotWidth );
  t1->SetBranchAddress("clusterstripChargeTotCharge",  &clusterstripChargeTotCharge );
  t1->SetBranchAddress("clusterstripChargeLocalTrackPhi",  &clusterstripChargeLocalTrackPhi );
  t1->SetBranchAddress("clusterstripChargeGlobalTrackPhi",  &clusterstripChargeGlobalTrackPhi );
  t1->SetBranchAddress("clusterstripChargeLocalTrackTheta",  &clusterstripChargeLocalTrackTheta );
  t1->SetBranchAddress("clusterstripChargeGlobalTrackTheta",  &clusterstripChargeGlobalTrackTheta );
  t1->SetBranchAddress("clusterstripChargeDetid",  &clusterstripChargeDetid );
  t1->SetBranchAddress("clusterstripChargetrackPt",  &clusterstripChargetrackPt );
  //  t1->SetBranchAddress("clusterstripChargeBdotY",  &clusterstripChargeBdotY );

  uint32_t evCount=0;
   
  cout << "in here a" << endl;
  Int_t nentries = (Int_t)t1->GetEntries();
  cout << "entries " << nentries << endl;
  cout << "in here b" << endl;

  ///fill variables from tree 1
  for (Int_t e=0; e<nentries; e++) 
    {
      t1->GetEntry(e);
          
      //per cluster
      uint32_t up = partition->size();
      for(uint32_t k=0; k<up;k++)
	{
	  if( partition->at(k) == sdId )
	    {
	      if(clusterlayerwheel->at(k) == layer)
		{
		  subpartition.push_back(partition->at(k));
		  subclustercharge.push_back(clustercharge->at(k));
		  subclusterwidth.push_back(clusterwidth->at(k));
		  subclusterlayerwheel.push_back(clusterlayerwheel->at(k));
		  subtsoslocalphi.push_back(tsoslocalphi->at(k));
		  subtsoslocaltheta.push_back(tsoslocaltheta->at(k));
		  subtsosBdotY.push_back(tsosBdotY->at(k));
		}
	    }
	}
      //per strip
      uint32_t upStrip = clusterstripChargeSubdetid->size();
      for(uint32_t k=0; k<upStrip ; k++)
	{
	  if( clusterstripChargeSubdetid->at(k) == sdId)
	    {
	      if(clusterstripChargeLayerwheel->at(k)== layer)
		{
		  subclusterstripChargeSubdetid.push_back(clusterstripChargeSubdetid->at(k));
		  subclusterstripCharge.push_back(clusterstripCharge->at(k));
		  subclusterstripChargeLayerwheel.push_back(clusterstripChargeLayerwheel->at(k));
		  subclusterstripChargeLocalTrackPhi.push_back(clusterstripChargeLocalTrackPhi->at(k));
		  subclusterstripChargeLocalTrackTheta.push_back(clusterstripChargeLocalTrackTheta->at(k));
		  subclusterstripChargeGlobalTrackTheta.push_back(clusterstripChargeGlobalTrackTheta->at(k));
		  subclusterstripChargeTotWidth.push_back(clusterstripChargeTotWidth->at(k));
		  subclusterstripChargeStripNr.push_back(clusterstripChargeStripNr->at(k));
		  subclusterstripChargeTotCharge.push_back(clusterstripChargeTotCharge->at(k));

		}
	    }
	}
    }

  TH1F* chargeForAllWidthsData = new TH1F("chargeForAllWidths", "chargeForAllWidths" , 100, 0, 1000 );
  TH1F* clusterAllWidthsData = new TH1F("clusterAllWidths", "chargeForAllWidths" , 20, 0, 20 );
  TH1F* clusterAllWidthsfinal = new TH1F("clusterAllWidthsfinal", "chargeForAllWidthsfinal" , 20, 0, 20 );

  for(uint32_t m = 0; m<subclustercharge.size(); m++)
    {
      chargeForAllWidthsData->Fill(subclustercharge.at(m));
      clusterAllWidthsData->Fill(subclusterwidth.at(m));
      int factor = subtsosBdotY.at(m) > 0 ? 1 : -1;
      if(factor*tan(subtsoslocaltheta.at(m))> 0.78 && subclusterwidth.at(m)> 3){
	clusterAllWidthsfinal->Fill(subclusterwidth.at(m));

      }
    }
  clusterAllWidthsfinal->SetMarkerStyle(kFullCircle);
  TCanvas c0("chargeForAllWidths","chargeForAllWidths");
  clusterAllWidthsfinal->DrawNormalized("P");
  c0.SaveAs("plot2.eps");

  vector<float> clusterVector;
  TProfile* leadingStripProfileData = new TProfile("leadingcluschargeData", "leadingcluschargeData", 550, 0, 550, 0, 300);
  TH1D* leadingStripSumData = new TH1D("leadingclusSumData", "leadingclusSumData", 800, 0, 800);
  float lstripNr = 0; 
  float lstripCh = 0; 
  float lClusWidth = 0 ; 
  float lClusCharge = 0 ; 
  uint32_t stripCounter = 0; 
  bool clusterEnd = true;  
  vector<double> stripChargeSum;
  stripChargeSum.resize(512,0);

  TProfile* clusterStripCahrgeAsFceTanThetaData = new TProfile("clusterStripCahrgeAsFceTanThetaData", "clusterStripCahrgeAsFceTanThetaData" , 50, -4, 4, 0, 300 );

  TProfile* clusterShapeData = new TProfile("clusterShapeData", "clusterShapeData", 40, -20, 20, 0, 300);
  TProfile* clusterShapeMT8Data = new TProfile("clusterShapeMT8Data", "clusterShapeMT8Data", 40, -20, 20, 0, 300);
  TProfile* clusterShapeLT7Data = new TProfile("clusterShapeML7Data", "clusterShapeLT7Data", 40, -20, 20, 0, 300);
  TProfile* clusterShape2Data = new TProfile("clusterShape2Data", "clusterShape2Data", 40, -20, 20, 0, 300);
  TProfile* clusterShape3Data = new TProfile("clusterShape3Data", "clusterShape3Data", 40, -20, 20, 0, 300);
  TProfile* clusterShape5Data = new TProfile("clusterShape5Data", "clusterShape5Data", 40, -20, 20, 0, 300);

  TProfile* clusterShapeDataPositive = new TProfile("clusterShapeDataPositive", "clusterShapeDataPositive", 40, 0, 20, 0, 300, "S");
  TProfile* clusterShapeMT8DataPositive = new TProfile("clusterShapeMT8DataPositive", "clusterShapeMT8DataPositive", 40, 0, 40, 0, 300, "S");
  TProfile* clusterShapeLT7DataPositive = new TProfile("clusterShapeML7DataPositive", "clusterShapeLT7DataPositive", 40, 0, 40, 0, 300, "S");
  TProfile* clusterShape2DataPositive = new TProfile("clusterShape2DataPositive", "clusterShape2DataPositive", 40, 0, 40, 0, 300, "S");
  TProfile* clusterShape3DataPositive = new TProfile("clusterShape3DataPositive", "clusterShape3DataPositive", 40, 0, 40, 0, 300, "S");

  TProfile* clusterShapeDataNegative = new TProfile("clusterShapeDataNegative", "clusterShapeDataNegative", 40, 0, 40, 0, 300, "S");
  TProfile* clusterShapeMT8DataNegative = new TProfile("clusterShapeMT8DataNegative", "clusterShapeMT8DataNegative", 40, 0, 40, 0, 300, "S");
  TProfile* clusterShapeLT7DataNegative = new TProfile("clusterShapeML7DataNegative", "clusterShapeLT7DataNegative", 40, 0, 40, 0, 300, "S");
  TProfile* clusterShape2DataNegative = new TProfile("clusterShape2DataNegative", "clusterShape2DataNegative", 40, 0, 40, 0, 300, "S");
  TProfile* clusterShape3DataNegative = new TProfile("clusterShape3DataNegative", "clusterShape3DataNegative", 40, 0, 40, 0, 300, "S");
  TH1F* clusterChargeWOSaturationData = new TH1F("clusterChargeWOSaturationData", "clusterChargeWOSaturationData", 100, 0, 1000 );
  TH1F* clusterWidthWOSaturationData = new TH1F("clusterWidthWOSaturationData", "clusterWidthWOSaturationData", 15, 0, 15 );
  TH1F* clusterChargePerStripWOSaturationData = new TH1F("clusterChargePerStripWOSaturationData", "clusterChargePerStripWOSaturationData", 300, 0, 300 );
  
  cout<<"in line 254"<<endl;
for(uint32_t m = 0; m<subclusterstripChargeLayerwheel.size(); m++)
    {
      if(subclusterstripChargeLayerwheel.at(m) == 3)
	{
	  //	  cout<<"in line 259 "<<endl;
	  if(lClusWidth==0 || clusterEnd==true )
	    {
	      cout<<"in line 262"<<endl;
	      lstripNr = subclusterstripChargeStripNr.at(m);
	      lstripCh = subclusterstripCharge.at(m);
	      lClusWidth = subclusterstripChargeTotWidth.at(m);
	      lClusCharge = subclusterstripChargeTotCharge.at(m);
	      cout << "cluster inital charge " << lstripCh << " for cluster width " << lClusWidth << "  m " << m << endl;
	      stripCounter = 1; 
	      clusterVector.clear();
	    }
	  if(stripCounter <= lClusWidth)
	    {
	      cout<<"in line 273"<<endl;
	      clusterEnd = false;
	      if(lstripCh<subclusterstripCharge.at(m))
		{
		  lstripNr = subclusterstripChargeStripNr.at(m);
		  lstripCh = subclusterstripCharge.at(m);
		 cout << "cluster charge changed to " << lstripCh << " m " << m << endl;
		}
	      clusterVector.push_back(subclusterstripCharge.at(m)) ;
	    
	  if(stripCounter == lClusWidth)
	    {
	      cout << "cluster charge final into profile " << lstripCh << " m " << m << endl;
	      leadingStripProfileData->Fill(lstripNr, lstripCh);
	      clusterStripCahrgeAsFceTanThetaData->Fill(tan(subclusterstripChargeLocalTrackTheta.at(m)), lstripCh);
	      double prev = stripChargeSum.at(lstripNr);
	      stripChargeSum.at(lstripNr) = prev+lstripCh;
	      float maxChrg = 0;
	      int32_t maxChrgCtr = -1;
	      for(uint32_t ch=0; ch<clusterVector.size();ch++)
		{
		  if(maxChrg < clusterVector.at(ch))
		    {
		      maxChrg = clusterVector.at(ch);
		      maxChrgCtr = ch;
		    }
		}
	      if(maxChrg < 253)
		{
		  clusterWidthWOSaturationData->Fill(lClusWidth);
		  clusterChargeWOSaturationData->Fill(lClusCharge);
		}
	      for(uint32_t ch=0; ch<clusterVector.size();ch++)
		{
		  if(maxChrg < 253)
		    clusterChargePerStripWOSaturationData->Fill(clusterVector.at(ch));
		  int32_t positionValue = ch-(maxChrgCtr);
		  cout << "charge position " << positionValue << " chareg " << clusterVector.at(ch) << endl;
		  clusterShapeData->Fill( positionValue ,clusterVector.at(ch) );
		  if( subclusterstripChargeLocalTrackTheta.at(m)>= 0)
		    clusterShapeDataPositive->Fill( ch+1 ,clusterVector.at(ch) );
		  if( subclusterstripChargeLocalTrackTheta.at(m)< 0)
		    clusterShapeDataNegative->Fill( ch+1 ,clusterVector.at(ch) );

		  if(lClusWidth>3)
		    {
		      clusterShape2Data->Fill( positionValue ,clusterVector.at(ch) );
		      if( subclusterstripChargeLocalTrackTheta.at(m)>= 0)
			clusterShape2DataPositive->Fill( ch+1 ,clusterVector.at(ch) );
		      if( subclusterstripChargeLocalTrackTheta.at(m)< 0)
			clusterShape2DataNegative->Fill( ch+1 ,clusterVector.at(ch) );
		    }
		  if(lClusWidth==3)
		    {
		      clusterShape3Data->Fill( positionValue ,clusterVector.at(ch) );
		      if( subclusterstripChargeLocalTrackTheta.at(m)>= 0)
			clusterShape3DataPositive->Fill( ch+1 ,clusterVector.at(ch) );
		      if( subclusterstripChargeLocalTrackTheta.at(m)< 0)
			clusterShape3DataNegative->Fill( ch+1 ,clusterVector.at(ch) );
		    }

		}
	      clusterVector.clear();
	      clusterEnd =true;
	    }
	  stripCounter++;
	    }
	}
    }

 clusterShape3DataPositive->SetMarkerStyle(kFullCircle); 
 clusterShape3DataNegative->SetMarkerStyle(kFullCircle); 
 clusterShape3DataPositive->SetMarkerColor(kPink); 
 clusterShape3DataPositive->SetMarkerColor(kBlue); 
 clusterShape3DataPositive->SetMaximum(1.5*  clusterShape3DataPositive->GetMaximum());
 TCanvas c22Positive("clusterShape3Positive","clusterShape3Positive");
 clusterShape3DataPositive->DrawNormalized("P"); 
 clusterShape3DataNegative->DrawNormalized("P same hist e"); 
 c22Positive.SaveAs("clustershape.pdf");

 clusterShape2DataPositive->SetMarkerStyle(kFullCircle);
 clusterShape2DataNegative->SetMarkerStyle(kFullCircle);
 clusterShape2DataPositive->SetMarkerColor(kPink);
 clusterShape2DataPositive->SetMarkerColor(kBlue);
 clusterShape2DataPositive->SetMaximum(1.5*  clusterShape2DataPositive->GetMaximum());
 TCanvas cPositive("clusterShape2Positive","clusterShape2Positive");
 clusterShape2DataPositive->DrawNormalized("P");
 clusterShape2DataNegative->DrawNormalized("P same hist e");
 cPositive.SaveAs("clustershapemorethan3.pdf");


 return 0;
}
