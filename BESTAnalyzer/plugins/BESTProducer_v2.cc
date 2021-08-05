// -*- C++ -*-
//
// Package:    VLQAnalysis/BESTAnalyzer
// Class:      BESTAnalyzer
//
/**\class BESTAnalyzer BESTAnalyzer.cc VLQAnalysis/BESTAnalyzer/plugins/BESTAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Reyer Band
//         Created:  Tue, 11 Feb 2020 19:27:20 GMT
//
//


// system include files
#include <memory>
#include <thread>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "PhysicsTools/CandUtils/interface/Thrust.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

// #include <fastjet/JetDefinition.hh>
// #include <fastjet/PseudoJet.hh>
// #include "fastjet/tools/Filter.hh"
// #include <fastjet/ClusterSequence.hh>
// #include <fastjet/ActiveAreaSpec.hh>
// #include <fastjet/ClusterSequenceArea.hh>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TCanvas.h"

//user made files
#include "VLQAnalysis/BESTAnalyzer/interface/BESTtoolbox.h"
#include "VLQAnalysis/BESTAnalyzer/interface/CacheHandler.h"
// #include "VLQAnalysis/BESTAnalyzer/src/CacheHandler.cc"
//???why include a .cc file?
#include "VLQAnalysis/BESTAnalyzer/interface/BESTEvaluation.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


//using reco::TrackCollection;

//class BESTAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
// class BESTAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns > {
class BESTProducer_v2 : public edm::stream::EDProducer<> {
    //???what is the difference of edm::stream::EDProducer<> than edm:::EDProducer 
public:
  //  explicit BESTAnalyzer(const edm::ParameterSet&, const CacheHandler*);
  explicit BESTProducer_v2(const edm::ParameterSet&);

  ~BESTProducer_v2();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  //  static std::unique_ptr<CacheHandler> initializeGlobalCache(const edm::ParameterSet& cfg);
  //  static void globalEndJob(const CacheHandler* cache) {}


private:
  // virtual void beginJob() override;
  // virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  // virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  // virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  // virtual void endJob() override;
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;


  // ----------member data ---------------------------
  //  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  //input variables from run config file
  std::string inputJetColl_;
  std::string GT_; //input global tag will be used to decide which year of MC/data is being input

  // Tree variables
  // TTree *jetTree;
  // TH1F *GenWeightTotal;
  // TH1F *Cutflow;
  std::map<std::string, float> treeVars;
  std::vector<std::string> listOfVars;
  std::map<std::string, std::vector<float> > treeVecVars;
  std::map<std::string, std::vector<int> > intVecVars;
  std::vector<std::string> listOfVecVars;
  std::vector<std::string> listOfIntVecVars;
  //Tokens
  edm::EDGetTokenT<std::vector<pat::Jet> > ak8JetsToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genPartToken_;
  edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate> > secVerticesToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken_;
  edm::EDGetTokenT<edm::TriggerResults> trigResultsToken_;
  edm::EDGetTokenT<bool> BadChCandFilterToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken_;
  edm::EDGetTokenT<LHEEventProduct> lheEventProductToken_;
  edm::EDGetTokenT<LHERunInfoProduct> lheRunInfoProductToken_;

  BESTEvaluation* BEST_;
  const CacheHandler* cache_;
  std::vector<std::string> listOfBESTVars_;
  std::string name_;
  edm::FileInPath path_;
  edm::FileInPath means_;
  bool isMC_;
  bool isSignal_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
//BESTAnalyzer::BESTAnalyzer(const edm::ParameterSet& iConfig, const CacheHandler* cache):
BESTProducer_v2::BESTProducer_v2(const edm::ParameterSet& iConfig):
  inputJetColl_ (iConfig.getParameter<std::string>("inputJetColl")),
  //  GT_ (iConfig.getParameter<std::string>("GT")),
  name_ (iConfig.getParameter<std::string>("name")), //BESTGraph
  path_ (iConfig.getParameter<edm::FileInPath>("path")), //VLQAnalysis/BESTAnalyzer/test/constantgraph.pb
  means_ (iConfig.getParameter<edm::FileInPath>("means")),//ScalerParameters.txt
  isMC_ (iConfig.getParameter<bool>("isMC")),
  isSignal_ (iConfig.getParameter<bool>("isSignal"))
{
  cache_ = new CacheHandler(path_); // Class to create an instance of a tensorflow session that loads the trained BEST network from a .pb file as a graph
  BEST_ = new BESTEvaluation(cache_); //Loads the BEST Neural Network using the tensforflow interface
  BEST_->configure(iConfig);



  //now do what ever initialization is needed
  //???not sure why the output file has a dir "run"?
  // edm::Service<TFileService> fs;
  // jetTree = fs->make<TTree>("jetTree","jetTree");
  // GenWeightTotal = fs->make<TH1F>("GenWeightTotal", "GenWeightTotal", 1, 0.5, 1.5);
  // Cutflow = fs->make<TH1F>("Cutflow", "Cutflow", 7, 0, 6);

  //Store the BEST variables for each jet
  listOfVars.push_back("nJets");//std::vector<std::string> listOfVars
  listOfVecVars.push_back("jetAK8_phi"); //std::vector<std::string>
  listOfVecVars.push_back("jetAK8_eta");
  listOfVecVars.push_back("jetAK8_pt");
  listOfVecVars.push_back("jetAK8_mass");
  listOfVecVars.push_back("jetAK8_SoftDropMass");
  listOfVecVars.push_back("nSecondaryVertices");
  listOfVecVars.push_back("bDisc");
  listOfVecVars.push_back("bDisc1");
  listOfVecVars.push_back("bDisc2");
  listOfVecVars.push_back("jetAK8_Tau4");
  listOfVecVars.push_back("jetAK8_Tau3");
  listOfVecVars.push_back("jetAK8_Tau2");
  listOfVecVars.push_back("jetAK8_Tau1");
  listOfVecVars.push_back("jetAK8_Tau32");
  listOfVecVars.push_back("jetAK8_Tau21");
  listOfVecVars.push_back("FoxWolfH1_Higgs");
  listOfVecVars.push_back("FoxWolfH2_Higgs");
  listOfVecVars.push_back("FoxWolfH3_Higgs");
  listOfVecVars.push_back("FoxWolfH4_Higgs");

  listOfVecVars.push_back("FoxWolfH1_Top");
  listOfVecVars.push_back("FoxWolfH2_Top");
  listOfVecVars.push_back("FoxWolfH3_Top");
  listOfVecVars.push_back("FoxWolfH4_Top");

  listOfVecVars.push_back("FoxWolfH1_W");
  listOfVecVars.push_back("FoxWolfH2_W");
  listOfVecVars.push_back("FoxWolfH3_W");
  listOfVecVars.push_back("FoxWolfH4_W");

  listOfVecVars.push_back("FoxWolfH1_Z");
  listOfVecVars.push_back("FoxWolfH2_Z");
  listOfVecVars.push_back("FoxWolfH3_Z");
  listOfVecVars.push_back("FoxWolfH4_Z");

  listOfVecVars.push_back("isotropy_Higgs");
  listOfVecVars.push_back("sphericity_Higgs");
  listOfVecVars.push_back("aplanarity_Higgs");
  listOfVecVars.push_back("thrust_Higgs");

  listOfVecVars.push_back("isotropy_Top");
  listOfVecVars.push_back("sphericity_Top");
  listOfVecVars.push_back("aplanarity_Top");
  listOfVecVars.push_back("thrust_Top");

  listOfVecVars.push_back("isotropy_W");
  listOfVecVars.push_back("sphericity_W");
  listOfVecVars.push_back("aplanarity_W");
  listOfVecVars.push_back("thrust_W");

  listOfVecVars.push_back("isotropy_Z");
  listOfVecVars.push_back("sphericity_Z");
  listOfVecVars.push_back("aplanarity_Z");
  listOfVecVars.push_back("thrust_Z");

  listOfVecVars.push_back("asymmetry_Higgs");
  listOfVecVars.push_back("asymmetry_Top");
  listOfVecVars.push_back("asymmetry_W");
  listOfVecVars.push_back("asymmetry_Z");

  listOfVecVars.push_back("nSubjets_Higgs");
  listOfVecVars.push_back("nSubjets_Top");
  listOfVecVars.push_back("nSubjets_W");
  listOfVecVars.push_back("nSubjets_Z");

  listOfVecVars.push_back("subjet12_mass_Higgs");
  listOfVecVars.push_back("subjet23_mass_Higgs");
  listOfVecVars.push_back("subjet13_mass_Higgs");
  listOfVecVars.push_back("subjet1234_mass_Higgs");

  listOfVecVars.push_back("subjet12_mass_Top");
  listOfVecVars.push_back("subjet23_mass_Top");
  listOfVecVars.push_back("subjet13_mass_Top");
  listOfVecVars.push_back("subjet1234_mass_Top");

  listOfVecVars.push_back("subjet12_mass_W");
  listOfVecVars.push_back("subjet23_mass_W");
  listOfVecVars.push_back("subjet13_mass_W");
  listOfVecVars.push_back("subjet1234_mass_W");

  listOfVecVars.push_back("subjet12_mass_Z");
  listOfVecVars.push_back("subjet23_mass_Z");
  listOfVecVars.push_back("subjet13_mass_Z");
  listOfVecVars.push_back("subjet1234_mass_Z");

  listOfVecVars.push_back("subjet12_CosTheta_Higgs");
  listOfVecVars.push_back("subjet23_CosTheta_Higgs");
  listOfVecVars.push_back("subjet13_CosTheta_Higgs");
  listOfVecVars.push_back("subjet1234_CosTheta_Higgs");

  listOfVecVars.push_back("subjet12_CosTheta_Top");
  listOfVecVars.push_back("subjet23_CosTheta_Top");
  listOfVecVars.push_back("subjet13_CosTheta_Top");
  listOfVecVars.push_back("subjet1234_CosTheta_Top");

  listOfVecVars.push_back("subjet12_CosTheta_W");
  listOfVecVars.push_back("subjet23_CosTheta_W");
  listOfVecVars.push_back("subjet13_CosTheta_W");
  listOfVecVars.push_back("subjet1234_CosTheta_W");

  listOfVecVars.push_back("subjet12_CosTheta_Z");
  listOfVecVars.push_back("subjet23_CosTheta_Z");
  listOfVecVars.push_back("subjet13_CosTheta_Z");
  listOfVecVars.push_back("subjet1234_CosTheta_Z");

  listOfVecVars.push_back("subjet12_DeltaCosTheta_Higgs");
  listOfVecVars.push_back("subjet13_DeltaCosTheta_Higgs");
  listOfVecVars.push_back("subjet23_DeltaCosTheta_Higgs");

  listOfVecVars.push_back("subjet12_DeltaCosTheta_Top");
  listOfVecVars.push_back("subjet13_DeltaCosTheta_Top");
  listOfVecVars.push_back("subjet23_DeltaCosTheta_Top");

  listOfVecVars.push_back("subjet12_DeltaCosTheta_W");
  listOfVecVars.push_back("subjet13_DeltaCosTheta_W");
  listOfVecVars.push_back("subjet23_DeltaCosTheta_W");

  listOfVecVars.push_back("subjet12_DeltaCosTheta_Z");
  listOfVecVars.push_back("subjet13_DeltaCosTheta_Z");
  listOfVecVars.push_back("subjet23_DeltaCosTheta_Z");

  //NN output
  listOfIntVecVars.push_back("BESTDecision");
  listOfVecVars.push_back("NNOutputs0");
  listOfVecVars.push_back("NNOutputs1");
  listOfVecVars.push_back("NNOutputs2");
  listOfVecVars.push_back("NNOutputs3");
  listOfVecVars.push_back("NNOutputs4");
  listOfVecVars.push_back("NNOutputs5");


  listOfVars.push_back("HT");
  listOfVars.push_back("HT__JEC_Up");
  listOfVars.push_back("HT__JEC_Dn");
  listOfVars.push_back("HT__JER_Up");
  listOfVars.push_back("HT__JER_Dn");

  //Trigger info For 2016
  listOfVars.push_back("PFHT800");
  listOfVars.push_back("PFHT900");
  listOfVars.push_back("PFJet400");
  //Trigger for 2017 and 2018
  listOfVars.push_back("PFHT1050");
  listOfVars.push_back("PFJet450");

  //MC Info
  listOfVars.push_back("EvtWeight");
  listOfVars.push_back("PileupWeight");
  listOfVars.push_back("PileupWeightUp");
  listOfVars.push_back("PileupWeightDn");
  listOfVars.push_back("PDFWeightUp");
  listOfVars.push_back("PDFWeightDn");
  listOfVars.push_back("Q2WeightUp");
  listOfVars.push_back("Q2WeightDn");

  //Signal MC Info
  listOfVars.push_back("VLQDecayMode");
  listOfIntVecVars.push_back("JetGenID");
  // for (unsigned i = 0; i < listOfVars.size(); i++){
    // treeVars[ listOfVars[i] ] = -999.99;
    // jetTree->Branch( (listOfVars[i]).c_str() , &(treeVars[ listOfVars[i] ]), (listOfVars[i]+"/F").c_str() );
  // }
//
  // for (unsigned i = 0; i < listOfVecVars.size(); i++){
    // jetTree->Branch( (listOfVecVars[i]).c_str() , &(treeVecVars[ listOfVecVars[i] ]) );
  // }
//
  // for (unsigned i = 0; i < listOfIntVecVars.size(); i++){
    // jetTree->Branch( (listOfIntVecVars[i]).c_str() , &(intVecVars[ listOfIntVecVars[i] ]) );
  // }

  listOfBESTVars_ = {"jetAK8_pt", "jetAK8_mass", "jetAK8_SoftDropMass", "nSecondaryVertices", "bDisc", "bDisc1", "bDisc2", "jetAK8_Tau4", "jetAK8_Tau3", "jetAK8_Tau2", "jetAK8_Tau1", "jetAK8_Tau32", "jetAK8_Tau21", "FoxWolfH1_Higgs", "FoxWolfH2_Higgs", "FoxWolfH3_Higgs", "FoxWolfH4_Higgs", "FoxWolfH1_Top", "FoxWolfH2_Top", "FoxWolfH3_Top", "FoxWolfH4_Top", "FoxWolfH1_W", "FoxWolfH2_W", "FoxWolfH3_W", "FoxWolfH4_W", "FoxWolfH1_Z", "FoxWolfH2_Z", "FoxWolfH3_Z", "FoxWolfH4_Z", "isotropy_Higgs", "sphericity_Higgs", "aplanarity_Higgs", "thrust_Higgs", "sphericity_Top", "aplanarity_Top", "thrust_Top", "sphericity_W", "aplanarity_W", "thrust_W", "sphericity_Z", "aplanarity_Z", "thrust_Z", "nSubjets_Higgs", "nSubjets_Top", "nSubjets_W", "nSubjets_Z", "subjet12_mass_Higgs", "subjet23_mass_Higgs", "subjet13_mass_Higgs", "subjet1234_mass_Higgs", "subjet12_mass_Top", "subjet23_mass_Top", "subjet13_mass_Top", "subjet1234_mass_Top", "subjet12_mass_W", "subjet23_mass_W", "subjet13_mass_W", "subjet1234_mass_W", "subjet12_mass_Z", "subjet23_mass_Z", "subjet13_mass_Z", "subjet1234_mass_Z", "subjet12_CosTheta_Higgs", "subjet23_CosTheta_Higgs", "subjet13_CosTheta_Higgs", "subjet1234_CosTheta_Higgs", "subjet12_CosTheta_Top", "subjet23_CosTheta_Top", "subjet13_CosTheta_Top", "subjet1234_CosTheta_Top", "subjet12_CosTheta_W", "subjet23_CosTheta_W", "subjet13_CosTheta_W", "subjet1234_CosTheta_W", "subjet12_CosTheta_Z", "subjet23_CosTheta_Z", "subjet13_CosTheta_Z", "subjet1234_CosTheta_Z", "subjet12_DeltaCosTheta_Higgs", "subjet13_DeltaCosTheta_Higgs", "subjet23_DeltaCosTheta_Higgs", "subjet12_DeltaCosTheta_Top", "subjet13_DeltaCosTheta_Top", "subjet23_DeltaCosTheta_Top", "subjet12_DeltaCosTheta_W", "subjet13_DeltaCosTheta_W", "subjet23_DeltaCosTheta_W", "subjet12_DeltaCosTheta_Z", "subjet13_DeltaCosTheta_Z","subjet23_DeltaCosTheta_Z", "asymmetry_Higgs", "asymmetry_Top", "asymmetry_W", "asymmetry_Z"};
  //------------------------------------------------------------------------------
  // Define input tags -----------------------------------------------------------
  //------------------------------------------------------------------------------

  // AK8 Jets
  edm::InputTag ak8JetsTag_;
  ak8JetsTag_ = edm::InputTag("slimmedJetsAK8", "", "PAT");
  ak8JetsToken_ = consumes<std::vector<pat::Jet> >(ak8JetsTag_);

  // Gen Particles

  if (isMC_){
    edm::InputTag genPartTag_;
    genPartTag_ = edm::InputTag("prunedGenParticles", "", "PAT");
    genPartToken_ = consumes<std::vector<reco::GenParticle> >(genPartTag_);
  }
  // Primary Vertices
  edm::InputTag verticesTag_;
  verticesTag_ = edm::InputTag("offlineSlimmedPrimaryVertices", "", "PAT");
  verticesToken_ = consumes<std::vector<reco::Vertex> >(verticesTag_);

  // Secondary Vertices
  edm::InputTag secVerticesTag_;
  secVerticesTag_ = edm::InputTag("slimmedSecondaryVertices", "", "PAT");
  secVerticesToken_ = consumes<std::vector<reco::VertexCompositePtrCandidate> >(secVerticesTag_);

  //Generator Event Info For Weights
  //Only called for MC
  if (isMC_){
    edm::InputTag genEvtInfoTag_;
    genEvtInfoTag_ = edm::InputTag("generator", "", "SIM");
    edm::InputTag lheRunInfoProductTag_;
    lheRunInfoProductTag_ = edm::InputTag("externalLHEProducer", "", "SIM");
    edm::InputTag lheEventProductTag_;
    lheEventProductTag_ = edm::InputTag("externalLHEProducer", "", "SIM");
    genEvtInfoToken_ = consumes<GenEventInfoProduct> (genEvtInfoTag_);
    lheRunInfoProductToken_ = consumes<LHERunInfoProduct, edm::InRun> (lheRunInfoProductTag_);
    lheEventProductToken_ = consumes<LHEEventProduct> (lheEventProductTag_);
  }

    produces<pat::JetCollection>();
    // produces < CollectionName > (label),When the class constructor is called, the framework should be instructed that  class will add something to the file

}

BESTProducer_v2::~BESTProducer_v2()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------

// void
// BESTProducer_v2::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
// {

// }
// BESTProducer_v2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
void BESTProducer_v2::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;


  Handle< std::vector<pat::Jet> > ak8JetsCollection;
  iEvent.getByToken(ak8JetsToken_, ak8JetsCollection);
  std::vector<pat::Jet> ak8Jets = *ak8JetsCollection.product();

  Handle< std::vector<reco::GenParticle> > genPartCollection;
  std::vector<reco::GenParticle> genPart;
  if (isMC_){
    iEvent.getByToken(genPartToken_, genPartCollection);
    genPart = *genPartCollection.product();
  }
  Handle< std::vector<reco::Vertex> > vertexCollection;
  iEvent.getByToken(verticesToken_, vertexCollection);
  std::vector<reco::Vertex> pVertices = *vertexCollection.product();

  Handle< std::vector<reco::VertexCompositePtrCandidate> > secVertexCollection;
  iEvent.getByToken(secVerticesToken_, secVertexCollection);
  std::vector<reco::VertexCompositePtrCandidate> secVertices = *secVertexCollection.product();

  Handle<GenEventInfoProduct> genEvtInfo;
  if(isMC_){
    iEvent.getByToken(genEvtInfoToken_, genEvtInfo);
  }

  //Get Generator Weights, Systematic Variations
  if(isMC_){
    float EventWeight = genEvtInfo->weight();
    // std::cout<<"EventWeight="<<EventWeight<<"\n";
    // GenWeightTotal->Fill(1, EventWeight);
    treeVars["EvtWeight"] = EventWeight;
  }
  //Cutflow: First stage will just equal to number of events
  // Cutflow->Fill(0);

  //Begin pre-selections
  auto outputs = std::make_unique<pat::JetCollection>();
  if (ak8Jets.size() > 3){


      // Cutflow->Fill(1);
      //Cutting on GeV > 400 for analysis, remember that the network is only trained on > 500!
      if (checkKinematicsOfJets(ak8Jets, 4) ){
          //in plugins/BESTtoolbox.cc
          // Cutflow->Fill(2);
          if (checkLengthOfSubJets(ak8Jets, 4) ){
              // Cutflow->Fill(3);

              treeVars["HT"] = ak8Jets[0].pt() + ak8Jets[1].pt() + ak8Jets[2].pt() + ak8Jets[3].pt();
              // std::cout<<"HT = "<<treeVars["HT"]<<"\n";


              //Fills map with basic kinematic variables
              //loop of jets begins
              for (int i = 0; i < 4; i++){
                  //???it seems it only uses leading 4 jet?
                  const pat::Jet& ijet = ak8Jets[i];
                  treeVecVars["jetAK8_phi"].push_back(ijet.phi());
                  treeVecVars["jetAK8_eta"].push_back(ijet.eta());





                  std::map<std::string, float> BESTmap;//so one map for ijet?
                  std::vector<float> BESTScores;

                  initBESTVars(BESTmap, listOfBESTVars_);
                  storeJetVariables(BESTmap, ijet, secVertices);

                  std::vector<reco::Candidate * > daughtersOfJet;
                  getJetDaughters(daughtersOfJet, ijet); //unzips the subjets and other daughters into one vector
                  if (daughtersOfJet.size() < 3) goto DontFill;
                  // Cutflow->Fill(4, 0.25); // 1/4 weight per jet

                  storeRestFrameVariables(BESTmap, daughtersOfJet, ijet, "Higgs", 125.);
                  if (BESTmap["nSubjets_Higgs"] < 3) goto DontFill; //Should be a cleaner way to check this before the first RestFrameVariables call

                  storeRestFrameVariables(BESTmap, daughtersOfJet, ijet, "Top", 172.5);
                  storeRestFrameVariables(BESTmap, daughtersOfJet, ijet, "W", 80.4);
                  storeRestFrameVariables(BESTmap, daughtersOfJet, ijet, "Z", 91.2);

                  // Cutflow->Fill(5, 0.25);

                  std::vector<float> BESTVars = orderBESTVars(BESTmap, listOfBESTVars_);


                  float HImage[31][31];
                  float TImage[31][31];
                  float WImage[31][31];
                  float ZImage[31][31];
                  prepareBoostedImage(ijet, daughtersOfJet, HImage, 125.);
                  prepareBoostedImage(ijet, daughtersOfJet, TImage, 172.5);
                  prepareBoostedImage(ijet, daughtersOfJet, WImage, 80.4);
                  prepareBoostedImage(ijet, daughtersOfJet, ZImage, 91.2);

                  for (auto it = BESTmap.cbegin(); it !=BESTmap.cend(); it++){
                      treeVecVars[it->first].push_back(it->second);
                  }
                  //Plug Jet values into network
                  BESTScores = BEST_->getPrediction(HImage,TImage,WImage,ZImage,BESTVars);
                  //Convert from a vector like (0,0,0,1,0) into an int with the decision
                  //Currently a float for dumb reasons
                  int decision = std::distance(BESTScores.begin(), std::max_element(BESTScores.begin(), BESTScores.end() ) );

                  treeVecVars["NNOutputs0"].push_back(BESTScores[0]);
                  treeVecVars["NNOutputs1"].push_back(BESTScores[1]);
                  treeVecVars["NNOutputs2"].push_back(BESTScores[2]);
                  treeVecVars["NNOutputs3"].push_back(BESTScores[3]);
                  treeVecVars["NNOutputs4"].push_back(BESTScores[4]);
                  treeVecVars["NNOutputs5"].push_back(BESTScores[5]);
                  intVecVars["BESTDecision"].push_back(decision);

                  if(isMC_){
                      intVecVars["JetGenID"].push_back(FindPDGid(ijet, genPart, isSignal_));
                  }


                //for output
                //for output jet
                // pat::Jet newJet( *ijet);
                pat::Jet newJet = ijet;
                std::cout<<"what's added to the newJet:"<<"\n";
                for (const auto &p : listOfVars){
                    newJet.addUserFloat("BEST_"+p, treeVars[ p ]);
                    std::cout<<p<<treeVars[p]<<" ";
                }
                std::cout<<"\n";
                //???how to add vector for each jet
                outputs->push_back(newJet);



              }//4 ijet loop
              // jetTree->Fill();
          }
      }
    iEvent.put(std::move(outputs));
                //???not  sure how output is saved into the event, which branch?
  }
  //-------------------------------------------------------------------------------
  // Clear and Reset all tree variables -------------------------------------------
  //-------------------------------------------------------------------------------
 DontFill: //Label to go-to when jets in loop fail a cut
  for (unsigned i = 0; i < listOfVars.size(); i++){
    treeVars[ listOfVars[i] ] = -999.99;
  }
  for (unsigned i = 0; i < listOfVecVars.size(); i++){
    treeVecVars[ listOfVecVars[i] ].clear();
  }
  for (unsigned i = 0; i < listOfIntVecVars.size(); i++){
    intVecVars[ listOfIntVecVars[i] ].clear();
  }


}


// ------------ method called once each job just before starting event loop  ------------
// void
// BESTProducer_v2::beginJob()
// {
// }

void
BESTProducer_v2::beginStream(edm::StreamID)
{
}
// ------------ method called once each job just after ending the event loop  ------------
// void
// BESTProducer_v2::endJob()
// {
// }

void
BESTProducer_v2::endStream()
{
}
// void
// BESTProducer_v2::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
// {
  /*
  edm::Handle<LHERunInfoProduct> lherun;
  typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;

  iRun.getByToken( lheRunInfoProductToken_, lherun );
  LHERunInfoProduct myLHERunInfoProduct = *(lherun.product());

  for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
    std::cout << iter->tag() << std::endl;
    std::vector<std::string> lines = iter->lines();
    for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
      std::cout << lines.at(iLine);
    }
  }
  */
// }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BESTProducer_v2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(BESTProducer_v2);
