#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <TLatex.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include <TString.h>
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooGaussModel.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TCut.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "RooHist.h"
#include "RooGenericPdf.h"
#include "RooTruthModel.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooProdPdf.h"
#include "RooEffProd.h"
#include <TRandom3.h>
#include "RooRandom.h"
#include "RooPoisson.h"
#include "RooGamma.h"
#include "RooBernstein.h"
#include "RooGaussian.h"
#include "iostream"
#include "fstream"
#include "RooFitResult.h"
#include "RooCBShape.h"
#include "RooConstVar.h"
#include "TNtupleD.h"

using namespace std;
using namespace RooFit;

void totalpdftoycon(){

  TRandom3 *rndGenerator = new TRandom3();
  int size=0,seed,unsucess=0;
  ofstream lfi;                                                                                                      
  TString title[10000];

  TNtupleD *nt = new TNtupleD("nt","","sig:sigpull:bd:bdpull:semi:bkg:taumean:tauerror:pull:taubd:taubderr:taubdpull");

  TNtupleD *ntlife = new TNtupleD("ntlife","","taumean:tauerror:pull");

  TNtupleD *nterrlife = new TNtupleD("nterrlife","","taumean:tauerror:pull");

  int signal,bdevnt,bkgevnt,semievnt;
  cout<<" Enter number of signal, Bd, semi,and background "<<endl;
  cin>>signal>>bdevnt>>semievnt>>bkgevnt;

  cout<<"size for fitting and seed "<<endl;
  cin>>size>>seed;
  for(int fc=0;fc<size;fc++){
    cout<<" \n acting on toy:"<<fc<<"\n\n"<<endl;
    
    RooRealVar* mass=new RooRealVar("mass","M_{#mu^{+}#mu^{-}}[ GeV/c^{2}]",4.9,5.9);
    RooRealVar* treco=new RooRealVar("treco","t_{reco}[ps]",1.2,12.5) ;
    RooRealVar* trecoe=new RooRealVar("trecoe","terr[ps]",0.013,0.26);
    
    //comb
    RooRealVar* B1=new RooRealVar("B1", "B_1comb", 0.516, 0.345, 0.687);
    RooFormulaVar* B2=new RooFormulaVar("B2", "B_2comb", "1.-@0", RooArgList(*B1));
    RooBernstein* mass_comb=new RooBernstein("mass_comb", "mass_comb", *mass, RooArgSet(*B1,*B2));
    RooRealVar* gB1=new RooRealVar("gB1", "B_1comb", 0.516);//, 0.345, 0.687);
    gB1->setConstant(true);
    RooFormulaVar* gB2=new RooFormulaVar("gB2", "B_2comb", "1.-@0", RooArgList(*gB1));
    RooBernstein* gmass_comb=new RooBernstein("gmass_comb", "mass_comb", *mass, RooArgSet(*gB1,*gB2));
    
    //Bsmass
    RooRealVar* meanBs=new RooRealVar("meanBs", "mean",5.3642,5.32, 5.40);
    RooRealVar *sigmacbbs=new RooRealVar("sigmacbbs","m1",0.0489);//,0.004,0.2);                                     
    sigmacbbs->setConstant(true);
    RooRealVar *tailcbbs=new RooRealVar("tailcbbs","",1.7617);//,0.1,10.0);                                            
    tailcbbs->setConstant(true);
    RooRealVar *powcbbs=new RooRealVar("powcbbs","",2.12);//,0,10);                                                 
    powcbbs->setConstant(true);
    RooCBShape *cbsbs=new RooCBShape("cbsbs","Signal Lineshape",*mass,*meanBs,*sigmacbbs,*tailcbbs,*powcbbs);

    RooRealVar* gmeanBs=new RooRealVar("gmeanBs", "mean",5.3642);//,5.35, 5.38);
    gmeanBs->setConstant(true);
    RooCBShape *gcbsbs=new RooCBShape("gcbsbs","Signal Lineshape",*mass,*gmeanBs,*sigmacbbs,*tailcbbs,*powcbbs);

    //Bdmass

    RooRealVar* meanBd=new RooRealVar("meanBd", "mean",5.2772);//,5.25, 5.30);
    meanBd->setConstant(true);
    RooRealVar* sigmaBd=new RooRealVar("sigmaBd","sigma1 of Gaussian",0.0379);//,0.0,0.4);
    sigmaBd->setConstant(true);
    RooGaussian* gausianBd=new RooGaussian("gausianBd","gaus1", *mass, *meanBd, *sigmaBd);

    RooRealVar *sigmacbbd=new RooRealVar("sigmacbbd","m1",0.0619);//,0.004,0.2);
    sigmacbbd->setConstant(true);
    RooRealVar *tailcbbd=new RooRealVar("tailcbbd","",1.4646);//,0.1,3.0);
    tailcbbd->setConstant(true);
    RooRealVar *powcbbd=new RooRealVar("powcbbd","",2.80);//,0,10);
    powcbbd->setConstant(true);
    RooCBShape *cbsbd=new RooCBShape("cbsbd","Signal Lineshape",*mass,*meanBd,*sigmacbbd,*tailcbbd,*powcbbd);
    RooRealVar* frac_gau_CB_bd=new RooRealVar("frac_gau_CB_bd","fraction", 0.556);//, 0, 1);
    frac_gau_CB_bd->setConstant(true);
    RooAddPdf* mass_pdf_bd=new RooAddPdf("mass_pdf_bd","gausian+CB shape", RooArgList(*gausianBd,*cbsbd), *frac_gau_CB_bd);
    
    //semi
    RooRealVar* meansemi=new RooRealVar("meansemi", "mean",4.868,4.7, 5.289);
    RooRealVar* sigmasemi=new RooRealVar("sigmasemi","sigma1 of Gaussian",0.155);//,0.0,1.0);
    RooGaussian* gausiansemi=new RooGaussian("gausiansemi","gaus1", *mass, *meansemi, *sigmasemi);
    sigmasemi->setConstant(true);
    RooRealVar* gmeansemi=new RooRealVar("gmeansemi", "mean",4.868);//,4.7, 5.289);
    gmeansemi->setConstant(true);
    RooGaussian* ggausiansemi=new RooGaussian("ggausiansemi","gaus1", *mass, *gmeansemi, *sigmasemi);
        

    RooRealVar* meant=new RooRealVar("meant", "mean", 0);
    RooRealVar* sigmat=new RooRealVar("sigmat","sigma",0.072);
    RooGaussModel* gm=new RooGaussModel("gm","gm", *treco,*meant,*sigmat);
  
    RooFormulaVar* effFor=new RooFormulaVar("effFor","-0.142-0.00483*treco+0.2056/(1+exp(-1.2822*treco))",RooArgSet(*treco));
    RooFormulaVar* effFor0=new RooFormulaVar("effFor0","0.0011+1/pow(treco,1.4821)-1.4374/(pow(treco,2)+0.51)",RooArgSet(*treco));
    
    RooTruthModel* CM =new RooTruthModel("CM","",*treco);

    //bd lifetime
    RooRealVar*  TauBd=new RooRealVar("TauBd","",1.45,0,100);

    RooDecay* Bkg_Ctaubd =new RooDecay("Bkg_Ctaubd","",*treco,*TauBd,*gm,RooDecay::SingleSided);
    RooFormulaVar* effForch1 = new RooFormulaVar("effForch1","8.95745e-07-4.95401e-08*treco-2.43181e-06/(1+exp(1.12513 *treco))", RooArgSet(*treco) );
    RooEffProd* CtEffBd =new RooEffProd("CtEffBd","",*Bkg_Ctaubd,*effForch1);

    RooRealVar*  gTauBd=new RooRealVar("gTauBd","",1.45);    
    gTauBd->setConstant(true);
    RooDecay* gBkg_Ctaubd =new RooDecay("gBkg_Ctaubd","",*treco,*gTauBd,*gm,RooDecay::SingleSided);
    RooEffProd* gCtEffBd =new RooEffProd("gCtEffBd","",*gBkg_Ctaubd,*effForch1);

    RooRealVar*  Tau=new RooRealVar("Tau","",1.70,0.00,10);//Bs lifetime
    RooDecay* sig_Ctau1 =new RooDecay("sig_Ctau1","",*treco,*Tau,*gm,RooDecay::SingleSided);
    RooRealVar*  gTau=new RooRealVar("gTau","",1.70);                              
    gTau->setConstant(true);
    RooDecay* gsig_Ctau1 =new RooDecay("gsig_Ctau1","",*treco,*gTau,*gm,RooDecay::SingleSided);

    //semi lifetime
    RooRealVar*  Tausemi=new RooRealVar("Tausemi","",1.338,0.0,10.0);
    RooDecay* semi_Ctau =new RooDecay("semi_Ctau","",*treco,*Tausemi,*gm,RooDecay::SingleSided);

    RooRealVar*  gTausemi=new RooRealVar("gTausemi","",1.338);
    gTausemi->setConstant(true);
    RooDecay* gsemi_Ctau =new RooDecay("gsemi_Ctau","",*treco,*gTausemi,*gm,RooDecay::SingleSided);

    // background lifetime                                                                                                    
    RooRealVar*  TauBkg2=new RooRealVar("TauBkg2","",3.1068,0,10);
    RooRealVar*  TauBkg3=new RooRealVar("TauBkg3","",0.2773,0,5);
    RooDecay* Bkg_Ctau3 =new RooDecay("Bkg_Ctau3","",*treco,*TauBkg3,*gm,RooDecay::SingleSided);
    RooDecay* Bkg_Ctau2 =new RooDecay("Bkg_Ctau2","",*treco,*TauBkg2,*gm,RooDecay::SingleSided);    
    RooRealVar* fg1=new RooRealVar("fg1","",0.588,0.0,1.0);
    RooAddPdf* bkg_Ctau=new RooAddPdf("bkg_Ctauf","",RooArgSet(*Bkg_Ctau2,*Bkg_Ctau3),*fg1);
    
    RooRealVar*  gTauBkg2=new RooRealVar("gTauBkg2","",3.1068);
    gTauBkg2->setConstant(true);
    RooRealVar*  gTauBkg3=new RooRealVar("gTauBkg3","",0.2773);
    gTauBkg3->setConstant(true);
    RooDecay* gBkg_Ctau3 =new RooDecay("gBkg_Ctau3","",*treco,*gTauBkg3,*gm,RooDecay::SingleSided);
    RooDecay* gBkg_Ctau2 =new RooDecay("gBkg_Ctau2","",*treco,*gTauBkg2,*gm,RooDecay::SingleSided);
    RooRealVar* gfg1=new RooRealVar("gfg1","",0.588);
    gfg1->setConstant(true);
    RooAddPdf* gbkg_Ctau=new RooAddPdf("gbkg_Ctauf","",RooArgSet(*gBkg_Ctau2,*gBkg_Ctau3),*gfg1);


    RooEffProd* CtEffSig =new RooEffProd("CtEffSig","",*sig_Ctau1,*effFor);
    RooEffProd* gCtEffSig =new RooEffProd("gCtEffSig","",*gsig_Ctau1,*effFor);
  
    RooRealVar* nbkg=new RooRealVar("nbkg","Background fraction",bkgevnt,20,550);
    RooRealVar* Nsig=new RooRealVar("Nsig","signal fraction",signal,0,150);
    RooRealVar* nbd1=new RooRealVar("nbd1","Bd",bdevnt,0,60);
    RooRealVar* nsemi=new RooRealVar("nsemi","semi",semievnt,0,500);
    
    RooRealVar* gnbkg=new RooRealVar("gnbkg","Background fraction",bkgevnt);
    RooRealVar* gNsig=new RooRealVar("gNsig","signal fraction",signal);
    RooRealVar* gnbd1=new RooRealVar("gnbd1","Bd",bdevnt);
    RooRealVar* gnsemi=new RooRealVar("gnsemi","semi",semievnt);

    int evnt=gnbkg->getVal()+gNsig->getVal()+gnbd1->getVal()+gnsemi->getVal();
    //==============
    RooProdPdf* BsPdf=new RooProdPdf ("BsPdf","",RooArgList(*cbsbs,*CtEffSig));
    RooProdPdf*  BdPdf=new RooProdPdf ("BdPdf","",RooArgList(*mass_pdf_bd,*CtEffBd) );
    RooProdPdf*  BkgPdf=new RooProdPdf ("BkgPdf","",RooArgList(*mass_comb,*bkg_Ctau) );
    RooProdPdf*  semiPdf=new RooProdPdf ("semiPdf","",RooArgList(*gausiansemi,*semi_Ctau) );
    

    RooProdPdf*  gBsPdf=new RooProdPdf ("gBsPdf","",RooArgList(*gcbsbs,*gCtEffSig));
    RooProdPdf*  gBdPdf=new RooProdPdf ("gBdPdf","",RooArgList(*mass_pdf_bd,*gCtEffBd));
    RooProdPdf*  gBkgPdf=new RooProdPdf ("gBkgPdf","",RooArgList(*gmass_comb,*gbkg_Ctau) );
    RooProdPdf*  gsemiPdf=new RooProdPdf ("gsemiPdf","",RooArgList(*ggausiansemi,*gsemi_Ctau) );


    cout<<"\n Experiment # "<<fc<<" with event "<<evnt<<"\n"<<endl;
    int seed1=(fc+1)*seed;
    
    RooAddPdf* TotPdf=new RooAddPdf("TotPdf"," Signal + Bkg Pdf",RooArgList(*BsPdf,*BdPdf,*BkgPdf,*semiPdf),RooArgList(*Nsig,*nbd1,*nbkg,*nsemi) );  
    RooAddPdf* gTotPdf=new RooAddPdf("gTotPdf"," Signal + Bkg Pdf",RooArgList(*gBsPdf,*gBdPdf,*gBkgPdf,*gsemiPdf),RooArgList(*gNsig,*gnbd1,*gnbkg,*gnsemi) );  
    RooGaussian* bdconstraint=new RooGaussian("bdconstraint","bd constraint",*nbd1,RooConst(bdevnt),RooConst(2)) ;
    RooGaussian* bdconstrainttau=new RooGaussian("bdconstrainttau","bd constraint",*TauBd,RooConst(1.45),RooConst(0.9)) ;
    RooGaussian* meanbsconstraint=new RooGaussian("meanbsconstraint","bd constraint",*meanBs,RooConst(5.36),RooConst(0.04)) ;
    RooGaussian* B1constraint=new RooGaussian("B1constraint", "constraint",*B1,RooConst(0.51),RooConst(0.2)) ;
    
    RooArgList constr;
    constr.add(*bdconstraint);
    constr.add(*meanbsconstraint);
    constr.add(*B1constraint);
    constr.add(*bdconstrainttau);
    
    RooProdPdf* constrpdf=new RooProdPdf("constrpdf","",constr);
    RooProdPdf* gmodel_con=new RooProdPdf("gmodel_con","model with constraint",RooArgSet(*gTotPdf,*constrpdf)) ;
    RooProdPdf* model_con=new RooProdPdf("model_con","model with constraint",RooArgSet(*TotPdf,*constrpdf)) ;

    RooDataSet* tmp_evt = bdconstraint->generate(RooArgSet(*nbd1),1);
    nbd1->setVal(tmp_evt->get(0)->getRealValue(nbd1->GetName()));
    delete tmp_evt;
    int bdevn=nbd1->getVal();
    
    RooDataSet* tmp_tau = bdconstrainttau->generate(RooArgSet(*TauBd),1);
    TauBd->setVal(tmp_tau->get(0)->getRealValue(TauBd->GetName()));
    delete tmp_tau;

    RooDataSet* tmp_evtB1 = B1constraint->generate(RooArgSet(*B1),1);
    B1->setVal(tmp_evtB1->get(0)->getRealValue(B1->GetName()));
    delete tmp_evtB1;

    RooDataSet* tmp_evtbs = meanbsconstraint->generate(RooArgSet(*meanBs),1);
    meanBs->setVal(tmp_evtbs->get(0)->getRealValue(meanBs->GetName()));
    delete tmp_evtbs;


    RooRandom::randomGenerator()->SetSeed(seed1);
    RooDataSet* data=model_con->generate(RooArgSet(*mass,*treco),evnt,Extended(true));           
    
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    RooMsgService::instance().setSilentMode(kTRUE);
    RooMsgService::instance().setStreamStatus(1,false);
    RooMsgService::instance().getStream(1).removeTopic(Integration) ;
    RooMsgService::instance().getStream(1).removeTopic(Minimization) ;
    RooMsgService::instance().getStream(1).removeTopic(Fitting) ;
    RooMsgService::instance().getStream(1).removeTopic(NumIntegration) ;
    RooMsgService::instance().getStream(1).removeTopic(Optimization) ;
    RooMsgService::instance().getStream(1).removeTopic(ObjectHandling) ;
    RooMsgService::instance().getStream(1).removeTopic(Eval) ;
    RooMsgService::instance().Print() ;    
    
    RooFitResult* fit3d=model_con->fitTo(*data,Extended(kTRUE),Save());//,ConditionalObservables(*trecoe));
    fit3d->Print("v");  
    double nsig=Nsig->getVal();
    double sigpl=(nsig-signal)/Nsig->getError();
    double nbd=nbd1->getVal();
    double Nsemi=nsemi->getVal();
    double bkgy=nbkg->getVal();
    double Ta=Tau->getVal();
    double taer=Tau->getError();
    double pul=(Ta-1.70)/taer;
    cout<<"pull "<<pul<<endl;
    int stat=fit3d->status();
    int cov=fit3d->covQual();
    double Taubd=TauBd->getVal();
    double Taubderr=TauBd->getError();

    if(stat==0&&cov==3){
      
      double result[12];
      result[0]=Nsig->getVal();
      result[1]=(Nsig->getVal()-signal)/Nsig->getError();
      result[2]=nbd1->getVal();
      result[3]=(nbd1->getVal()-bdevnt)/nbd1->getError();
      result[4]=nsemi->getVal();
      result[5]=nbkg->getVal();
      result[6]=Tau->getVal();
      result[7]=Tau->getError();
      result[8]=(Tau->getVal()-1.70)/Tau->getError();
      result[9]=TauBd->getVal();
      result[10]=TauBd->getError();
      result[11]=(TauBd->getVal()-1.60)/TauBd->getError();
      nt->Fill(result);
      if(Ta>3){   
      double study[3];
      study[0]=Tau->getVal();
      study[1]=Tau->getError();
      study[2]=(Tau->getVal()-1.70)/Tau->getError();
      /*
      RooPlot* mframe  = mass->frame(Title("B_{s} mass distribution"),Bins(50));                                             
      data->plotOn(mframe);                                                                                               
      TotPdf->plotOn(mframe);                                                                                                 
     TotPdf->plotOn(mframe, Components(*mass_comb), LineStyle(kDashed),LineColor(kGreen));                                 
     TotPdf->plotOn(mframe, Components(*mass_pdf_bd), LineStyle(kDashed),LineColor(kCyan));                             
     TotPdf->plotOn(mframe, Components(*cbsbs), LineStyle(kDashed),LineColor(kRed));                                
     TotPdf->plotOn(mframe, Components(*gausiansemi), LineStyle(kDashed),LineColor(kMagenta));                            
     mframe->Draw();                                                                                          
     TCanvas* canm=new TCanvas(Form("canm_%d",fc),"",600,600);                                                       
     mframe->Draw();
     TPaveText* paveText16mch0 = new TPaveText(0.62,0.70,0.80,0.88,"NDC");
     paveText16mch0->SetBorderSize(0.0);
     paveText16mch0->SetFillColor(kWhite);
     paveText16mch0->SetFillStyle(0);
     paveText16mch0->SetTextSize(0.02);
     paveText16mch0->AddText(Form("Nsig = %.2f #pm %.2f ", Nsig->getVal(), Nsig->getError()));
     paveText16mch0->AddText(Form("nbkg = %.2f #pm %.2f ", nbkg->getVal(), nbkg->getError()));
     paveText16mch0->AddText(Form("nsemi= %.2f #pm %.2f ", nsemi->getVal(), nsemi->getError()));
     paveText16mch0->AddText(Form("nbd1 = %.2f #pm %.2f ", nbd1->getVal(), nbd1->getError()));
     paveText16mch0->Draw();
     canm->SaveAs(Form("mass_2d_toy_%d.pdf",fc));                                                                     
     RooPlot* tframe  = treco->frame(Title("Decay time"),Bins(40));                                                 
     data->plotOn(tframe);                                                                                                 
     TotPdf->plotOn(tframe);                                                                                         
     TotPdf->plotOn(tframe, Components(*CtEffSig), LineStyle(kDashed),LineColor(kRed));                                   
     TotPdf->plotOn(tframe, Components(*bkg_Ctau), LineStyle(kDashed),LineColor(kGreen));                            
     TotPdf->plotOn(tframe, Components(*Bkg_Ctaubd), LineStyle(kDashed),LineColor(kCyan));                              
     TotPdf->plotOn(tframe, Components(*semi_Ctau), LineStyle(kDashed),LineColor(kMagenta));                         
     TCanvas* tmm=new TCanvas(Form("tmm_%d",fc),"",600,600);                                                          
     TPaveText* pav0 = new TPaveText(0.62,0.70,0.80,0.88,"NDC");
     pav0->SetBorderSize(0.0);
     pav0->SetFillColor(kWhite);
     pav0->SetFillStyle(0);
     pav0->SetTextSize(0.02);
     pav0->AddText(Form("#tau = %.3f #pm %.3f ps ", Tau.getVal(), Tau.getError()));
     pav0->AddText(Form("#tau_{semi} = %.3f #pm %.3f ps ", Tausemi.getVal(), Tausemi.getError()));
     pav0->AddText(Form("#tau_{B_{d}} = %.3f #pm %.3f ps ", TauBd.getVal(), TauBd.getError()));
     tframe->Draw();
     pav0->Draw();
     tmm->SetLogy();				
     tmm->SaveAs(Form("treco_2d_toy_%d.pdf",fc)); 
      */ntlife->Fill(study);
      }
      if(taer>1){
	double studyerr[3];
        studyerr[0]=Tau->getVal();
	studyerr[1]=Tau->getError();
        studyerr[2]=(Tau->getVal()-1.70)/Tau->getError();
	nterrlife->Fill(studyerr);
      }

    }
    delete data;
    /*if(pul>0.5&&pul<1){
      RooPlot* mframe  = mass->frame(Title("B_{s} mass distribution"),Bins(120));
      data->plotOn(mframe);
      TotPdf->plotOn(mframe);
  TotPdf->plotOn(mframe, Components(*mass_comb), LineStyle(kDashed),LineColor(kGreen));
  TotPdf->plotOn(mframe, Components(*mass_pdf_bd), LineStyle(kDashed),LineColor(kCyan));
  TotPdf->plotOn(mframe, Components(*mass_pdf_bs), LineStyle(kDashed),LineColor(kRed));
  TotPdf->plotOn(mframe, Components(*gausiansemi), LineStyle(kDashed),LineColor(kMagenta));

  TCanvas* canm=new TCanvas(Form("canm_%d",fc),"",600,600);
  mframe->Draw();
  canm->SaveAs(Form("mass_2d_toy_%d.pdf",fc));
  RooPlot* tframe  = treco->frame(Title("Decay time"),Bins(40));
  data->plotOn(tframe);
  TotPdf->plotOn(tframe);
  TotPdf->plotOn(tframe, Components(*CtEffSig), LineStyle(kDashed),LineColor(kRed));
  TotPdf->plotOn(tframe, Components(*bkg_Ctau), LineStyle(kDashed),LineColor(kGreen));
  TotPdf->plotOn(tframe, Components(*Bkg_Ctaubd), LineStyle(kDashed),LineColor(kCyan));
  TotPdf->plotOn(tframe, Components(*semi_Ctau), LineStyle(kDashed),LineColor(kMagenta));
  TCanvas* tmm=new TCanvas(Form("tmm_%d",fc),"",600,600);
  tframe->Draw();               
  tmm->SetLogy();                                                                                                                       
  tmm->SaveAs(Form("treco_2d_toy_%d.pdf",fc));
  }*/
  }
  TFile *f = new TFile(Form("2Dtoy_new_consgen_bsB1constraint_%d_%d.root",signal,bdevnt),"recreate");
  nt->Write();
  ntlife->Write();
  nterrlife->Write();

}
