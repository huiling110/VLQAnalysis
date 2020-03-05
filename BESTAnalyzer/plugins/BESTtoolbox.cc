//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// BESTtoolbox.cpp -------------------------------------------------------------------
//========================================================================================
// C++ file containing functions for use with CMS EDAnalyzer and EDProducer --------------
//////////////////////////////////////////////////////////////////////////////////////////

#include "VLQAnalysis/BESTAnalyzer/interface/BESTtoolbox.h"

//========================================================================================
// Calculate Legendre Polynomials --------------------------------------------------------
//----------------------------------------------------------------------------------------
// Simple Legendre polynomial function that can calculate up to order 4 ------------------
// Inputs: argument of the polynomial and order desired ----------------------------------
//----------------------------------------------------------------------------------------

float LegendreP(float x, int order){
   if (order == 0) return 1;
   else if (order == 1) return x;
   else if (order == 2) return 0.5*(3*x*x - 1);
   else if (order == 3) return 0.5*(5*x*x*x - 3*x);
   else if (order == 4) return 0.125*(35*x*x*x*x - 30*x*x + 3);
   else return 0;
}

//========================================================================================
// Calculate Fox Wolfram Moments ---------------------------------------------------------
//----------------------------------------------------------------------------------------
// This function calculates the Fox Wolfram moments for jet constituents -----------------
// in various rest frames. ---------------------------------------------------------------
// Inputs: particles (jet constiuents boosted to rest frame) and empty array that --------
//         that will store the FW moments ------------------------------------------------
//----------------------------------------------------------------------------------------

int FWMoments(std::vector<TLorentzVector> particles, double (&outputs)[5] ){

   // get number of particles to loop over
   int numParticles = particles.size();

   // get energy normalization for the FW moments
   float s = 0.0;
   for(int i = 0; i < numParticles; i++){
   	s += particles[i].E();
   }

   float H0 = 0.0;
   float H4 = 0.0;
   float H3 = 0.0;
   float H2 = 0.0;
   float H1 = 0.0;

   for (int i = 0; i < numParticles; i++){

   	for (int j = i; j < numParticles; j++){

                // calculate cos of jet constituent angles
   		float costh = ( particles[i].Px() * particles[j].Px() + particles[i].Py() * particles[j].Py()
                                   + particles[i].Pz() * particles[j].Pz() ) / ( particles[i].P() * particles[j].P() );
   		float w1 = particles[i].P();
   		float w2 = particles[j].P();

                // calculate legendre polynomials of jet constiteuent angles
   		float fw0 = LegendreP(costh, 0);
   		float fw1 = LegendreP(costh, 1);
   		float fw2 = LegendreP(costh, 2);
   		float fw3 = LegendreP(costh, 3);
   		float fw4 = LegendreP(costh, 4);

                // calculate the Fox Wolfram moments
   		H0 += w1 * w2 * fw0;
   		H1 += w1 * w2 * fw1;
   		H2 += w1 * w2 * fw2;
   		H3 += w1 * w2 * fw3;
   		H4 += w1 * w2 * fw4;

   	}
   }

   // Normalize the Fox Wolfram moments
   if (H0 == 0) H0 += 0.001;      // to prevent dividing by zero
   outputs[0] = (H0);
   outputs[1] = (H1 / H0);
   outputs[2] = (H2 / H0);
   outputs[3] = (H3 / H0);
   outputs[4] = (H4 / H0);

   return 0;
}

//========================================================================================
// Get All Jet Constituents --------------------------------------------------------------
//----------------------------------------------------------------------------------------
// This gets all the jet constituents (daughters) and stores them as a standard ----------
// vector --------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void getJetDaughters(std::vector<reco::Candidate * > &daughtersOfJet, const pat::Jet &jet){
    // First get all daughters for the first Soft Drop Subjet
    for (unsigned int i = 0; i < jet.daughter(0)->numberOfDaughters(); i++){
      if (jet.daughter(0)->daughter(i)->pt() < 0.5) continue;
      daughtersOfJet.push_back( (reco::Candidate *) jet.daughter(0)->daughter(i) );
    }
    // Get all daughters for the second Soft Drop Subjet
    for (unsigned int i = 0; i < jet.daughter(1)->numberOfDaughters(); i++){
      if (jet.daughter(1)->daughter(i)->pt() < 0.5) continue;
      daughtersOfJet.push_back( (reco::Candidate *) jet.daughter(1)->daughter(i));
    }
    // Get all daughters not included in Soft Drop
    for (unsigned int i = 2; i< jet.numberOfDaughters(); i++){
      if (jet.daughter(i)->pt() < 0.5) continue;
      daughtersOfJet.push_back( (reco::Candidate *) jet.daughter(i) );
    }
}

//========================================================================================
// Store Jet Variables -------------------------------------------------------------------
//----------------------------------------------------------------------------------------
// This takes various jet quantaties and stores them on the map used to fill -------------
// the jet tree --------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void storeJetVariables(std::map<std::string, float> &BESTVars, const pat::Jet &jet, const std::vector<reco::VertexCompositePtrCandidate> &secVertices)
{//This method stores the lab frame variables used in BEST to a map
  
  BESTVars["jetAK8_pt"] = jet.pt();
  BESTVars["jetAK8_mass"] = jet.mass();
  BESTVars["jetAK8_SoftDropMass"] = jet.userFloat("ak8PFJetsPuppiSoftDropMass");
  BESTVars["jetAK8_Tau4"] = jet.userFloat("NjettinessAK8Puppi:tau4");
  BESTVars["jetAK8_Tau3"] = jet.userFloat("NjettinessAK8Puppi:tau3");
  BESTVars["jetAK8_Tau2"] = jet.userFloat("NjettinessAK8Puppi:tau2");
  BESTVars["jetAK8_Tau1"] = jet.userFloat("NjettinessAK8Puppi:tau1");
  BESTVars["jetAK8_Tau21"] = jet.userFloat("NjettinessAK8Puppi:tau2") / jet.userFloat("NjettinessAK8Puppi:tau1");
  BESTVars["jetAK8_Tau32"] = jet.userFloat("NjettinessAK8Puppi:tau3") / jet.userFloat("NjettinessAK8Puppi:tau2");
  BESTVars["jetAK8_SoftDropMass"] = jet.userFloat("ak8PFJetsPuppiSoftDropMass");
  BESTVars["bDisc"] = jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb");
  const auto subjets = jet.subjets("SoftDropPuppi");
  BESTVars["bDisc1"] = subjets[0]->bDiscriminator("pfDeepCSVJetTags:probb") + subjets[0]->bDiscriminator("pfDeepCSVJetTags:probbb");
  BESTVars["bDisc2"] = subjets[1]->bDiscriminator("pfDeepCSVJetTags:probb") + subjets[1]->bDiscriminator("pfDeepCSVJetTags:probbb");

  int numMatched = 0;
  for(std::vector<reco::VertexCompositePtrCandidate>::const_iterator vertBegin = secVertices.begin(),
	vertEnd = secVertices.end(), ivert = vertBegin; ivert != vertEnd; ivert++){
    TLorentzVector vert(ivert->px(), ivert->py(), ivert->pz(), ivert->energy() );
    TLorentzVector ijet(jet.px(), jet.py(), jet.pz(), jet.energy());
    if(ijet.DeltaR(vert) < 0.8 ){
      numMatched++;
    }
  }
  BESTVars["nSecondaryVertices"] = numMatched;  

}

//========================================================================================
// Store Secondary Vertex Information ----------------------------------------------------
//----------------------------------------------------------------------------------------
// This takes various secondary vertex quantities and stores them on the map -------------
// used to fill the tree -----------------------------------------------------------------
//----------------------------------------------------------------------------------------

void storeSecVertexVariables(std::map<std::string, float> &BESVars,
                             std::map<std::string, std::vector<float> > &jetVecVars, TLorentzVector jet,
                             std::vector<reco::VertexCompositePtrCandidate> secVertices){

   int numMatched = 0; // counts number of secondary vertices
   for(std::vector<reco::VertexCompositePtrCandidate>::const_iterator vertBegin = secVertices.begin(),
              vertEnd = secVertices.end(), ivert = vertBegin; ivert != vertEnd; ivert++){
      TLorentzVector vert(ivert->px(), ivert->py(), ivert->pz(), ivert->energy() );
      // match vertices to jet
      if(jet.DeltaR(vert) < 0.8 ){
         numMatched++;
         // save secondary vertex info for the first three sec vertices
         jetVecVars["SV_pt"].push_back(ivert->pt() );
         jetVecVars["SV_eta"].push_back(ivert->eta() );
         jetVecVars["SV_phi"].push_back(ivert->phi() );
         jetVecVars["SV_mass"].push_back(ivert->mass() );
         jetVecVars["SV_nTracks"].push_back(ivert->numberOfDaughters() );
         jetVecVars["SV_chi2"].push_back(ivert->vertexChi2() );
         jetVecVars["SV_Ndof"].push_back(ivert->vertexNdof() );
      }
   }
   BESVars["nSecondaryVertices"] = numMatched;
}

//========================================================================================
// Store Rest Frame Variables ------------------------------------------------------------
//----------------------------------------------------------------------------------------
// This boosts an ak8 jet (and all of its constituents) into heavy object rest frame -----
// and then uses it to calculate FoxWolfram moments, Event Shape Variables, --------------
// and assymmetry variables --------------------------------------------------------------
//----------------------------------------------------------------------------------------

void storeRestFrameVariables(std::map<std::string, float> &BESTVars, std::vector<reco::Candidate *> daughtersOfJet,
			     const pat::Jet &jet, std::string frame, float mass){

    // get 4 vector for heavy object rest frame
    typedef reco::Candidate::PolarLorentzVector fourv;
    fourv thisJet = jet.polarP4();
    TLorentzVector thisJetLV(0.,0.,0.,0.);
    thisJetLV.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), mass );

    std::vector<TLorentzVector> particles;
    std::vector<math::XYZVector> particles2;
    std::vector<reco::LeafCandidate> particles3;
    std::vector<fastjet::PseudoJet> FJparticles;

    // 4 vectors to be filled with subjet additions
    TLorentzVector subjet12LV(0.,0.,0.,0.);
    TLorentzVector subjet13LV(0.,0.,0.,0.);
    TLorentzVector subjet23LV(0.,0.,0.,0.);
    TLorentzVector subjet1234LV(0.,0.,0.,0.);

    double sumPz = 0;
    double sumP = 0;

    // Boost to object rest frame
    for(unsigned int i = 0; i < daughtersOfJet.size(); i++){
        // Do not include low pT particles
        if (daughtersOfJet[i]->pt() < 0.5) continue;

        // Create 4 vector to boost to Higgs frame
        TLorentzVector thisParticleLV( daughtersOfJet[i]->px(), daughtersOfJet[i]->py(), daughtersOfJet[i]->pz(), daughtersOfJet[i]->energy() );

        // Boost to Higgs rest frame
        thisParticleLV.Boost( -thisJetLV.BoostVector() );


        // Now that PF candidates are stored, make the boost axis the Z-axis
        // Important for BES variables
        pboost( thisJetLV.Vect(), thisParticleLV.Vect(), thisParticleLV);

        particles.push_back( thisParticleLV );
        particles2.push_back( math::XYZVector( thisParticleLV.X(), thisParticleLV.Y(), thisParticleLV.Z() ));
        particles3.push_back( reco::LeafCandidate(+1, reco::Candidate::LorentzVector( thisParticleLV.X(), thisParticleLV.Y(),
                                                                                        thisParticleLV.Z(), thisParticleLV.T() ) ));
        FJparticles.push_back( fastjet::PseudoJet( thisParticleLV.X(), thisParticleLV.Y(), thisParticleLV.Z(), thisParticleLV.T() ) );

        // Sum rest frame momenta for asymmetry calculation, but only for pt > 10
	//Why?????
	//        if (daughtersOfJet[i]->pt() < 10) continue;
        sumPz += thisParticleLV.Pz();
        sumP += abs( thisParticleLV.P() );
    }

    // Fox Wolfram Moments
    double fwm[5] = { 0.0, 0.0 ,0.0 ,0.0,0.0};
    FWMoments( particles, fwm);
    BESTVars["FoxWolfH1_"+frame] = fwm[1];
    BESTVars["FoxWolfH2_"+frame] = fwm[2];
    BESTVars["FoxWolfH3_"+frame] = fwm[3];
    BESTVars["FoxWolfH4_"+frame] = fwm[4];

    // Event Shape Variables
    EventShapeVariables eventShapes( particles2 );
    Thrust thrustCalculator( particles3.begin(), particles3.end() );
    BESTVars["isotropy_"+frame]   = eventShapes.isotropy();
    BESTVars["sphericity_"+frame] = eventShapes.sphericity();
    BESTVars["aplanarity_"+frame] = eventShapes.aplanarity();
    BESTVars["thrust_"+frame]     = thrustCalculator.thrust();

    // Jet Asymmetry
    double asymmetry             = sumPz/sumP;
    BESTVars["asymmetry_"+frame] = asymmetry;

    // Recluster the jets in the heavy object rest frame
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
    fastjet::ClusterSequence cs(FJparticles, jet_def);
    //    std::vector<fastjet::PseudoJet> jetsFJ = sorted_by_pt(cs.inclusive_jets(20.0));
    //Changed to 0.0 here, the 20.0 cuts on pT relative to a meaningless axis
    //    std::vector<fastjet::PseudoJet> jetsFJ = sorted_by_pt(cs.inclusive_jets(0.0));
    std::vector<fastjet::PseudoJet> jetsFJ = sorted_by_E(cs.inclusive_jets(0.0));

    // Store reclustered jet info
    for(unsigned int i = 0; i < jetsFJ.size(); i++){
        // make a TLorentzVector for the current clustered subjet
        TLorentzVector iSubjetLV(jetsFJ[i].px(), jetsFJ[i].py(), jetsFJ[i].pz(), jetsFJ[i].e() );

        // get subjet four vector combinations
        switch(i){
        case 0:
            subjet12LV   = subjet12LV   + iSubjetLV;
            subjet13LV   = subjet13LV   + iSubjetLV;
            subjet1234LV = subjet1234LV + iSubjetLV;
            break;
        case 1:
	  BESTVars["subjet12_DeltaCosTheta_"+frame]   = subjet12LV.CosTheta() - iSubjetLV.CosTheta();
            subjet12LV   = subjet12LV   + iSubjetLV;
            subjet23LV   = subjet23LV   + iSubjetLV;
            subjet1234LV = subjet1234LV + iSubjetLV;
            break;
        case 2:
	  BESTVars["subjet13_DeltaCosTheta_"+frame]   = subjet13LV.CosTheta() - iSubjetLV.CosTheta();
            subjet13LV   = subjet13LV   + iSubjetLV;
	    BESTVars["subjet23_DeltaCosTheta_"+frame]   = subjet23LV.CosTheta() - iSubjetLV.CosTheta();
            subjet23LV   = subjet23LV   + iSubjetLV;
            subjet1234LV = subjet1234LV + iSubjetLV;
            break;
        case 3:
            subjet1234LV = subjet1234LV + iSubjetLV;
            break;
        }
    }

    // Store subjet mass combinations
    BESTVars["subjet12_mass_"+frame]   = subjet12LV.M();
    BESTVars["subjet13_mass_"+frame]   = subjet13LV.M();
    BESTVars["subjet23_mass_"+frame]   = subjet23LV.M();
    BESTVars["subjet1234_mass_"+frame] = subjet1234LV.M();
    BESTVars["subjet12_CosTheta_"+frame]   = subjet12LV.CosTheta();
    BESTVars["subjet13_CosTheta_"+frame]   = subjet13LV.CosTheta();
    BESTVars["subjet23_CosTheta_"+frame]   = subjet23LV.CosTheta();
    BESTVars["subjet1234_CosTheta_"+frame] = subjet1234LV.CosTheta();
    BESTVars["nSubjets_"+frame] = jetsFJ.size();
}

//========================================================================================
// Make boost axis the rest frame z axis -------------------------------------------------
//----------------------------------------------------------------------------------------
// Given jet constituent lab momentum, find momentum relative to beam direction pbeam ----
// plab = Particle 3-vector in Boost Frame -----------------------------------------------
// pbeam = Lab Jet 3-vector --------------------------------------------------------------
//----------------------------------------------------------------------------------------

void pboost( TVector3 pbeam, TVector3 plab, TLorentzVector &pboo ){

    double pl = plab.Dot(pbeam);
    pl *= double(1. / pbeam.Mag());

    // set x axis direction along pbeam x (0,0,1)
    TVector3 pbx;

    pbx.SetX(pbeam.Y());
    pbx.SetY(-pbeam.X());
    pbx.SetZ(0.0);

    pbx *= double(1. / pbx.Mag());

    // set y axis direction along -pbx x pbeam
    TVector3 pby;

    pby = -pbx.Cross(pbeam);
    pby *= double(1. / pby.Mag());

    pboo.SetX((plab.Dot(pbx)));
    pboo.SetY((plab.Dot(pby)));
    pboo.SetZ(pl);

}

void initBESTVars(std::map<std::string, float> &BESMap, std::vector<std::string> listOfBESTVars){
  //Takes an empty map and fills it in the correct order with the BEST inputs
  //This order is determined by the variable order presented in the network training
  for (unsigned int i=0; i < listOfBESTVars.size(); i++){
    BESMap[listOfBESTVars[i]] = -999;
  }
}

std::vector<float> orderBESTVars(std::map<std::string, float> BESMap, std::vector<std::string> listOfBESTVars){
  std::vector<float> orderedBESTVars;
  for (unsigned int i=0; i < listOfBESTVars.size(); i++){
    orderedBESTVars.push_back(BESMap[listOfBESTVars[i]]);
  }
  return orderedBESTVars;
}

void prepareBoostedImage(const pat::Jet &jet, const std::vector<reco::Candidate *> daughtersOfJet, float Image[31][31], const float mass){
  TLorentzVector thisJetLV(0.,0.,0.,0.);
  reco::Candidate::PolarLorentzVector thisJet = jet.polarP4();
  thisJetLV.SetPtEtaPhiM(thisJet.Pt(), thisJet.Eta(), thisJet.Phi(), mass );
  std::vector<TLorentzVector>* BoostedDaughters = new std::vector<TLorentzVector>;

  float sumE = 0;
  //Create a vector of candidates that is boosted, low pT filtered
  for(unsigned int i = 0; i < daughtersOfJet.size(); i++){
    if (daughtersOfJet[i]->pt() < 0.5) continue;
    TLorentzVector thisParticleLV( daughtersOfJet[i]->px(), daughtersOfJet[i]->py(), daughtersOfJet[i]->pz(), daughtersOfJet[i]->energy() );
    thisParticleLV.Boost( -thisJetLV.BoostVector() );
    sumE+= thisParticleLV.E();
    BoostedDaughters->push_back(thisParticleLV);
  }


  //Sort the new list by energy
  auto sortLambda = [] (const TLorentzVector& lv1, const TLorentzVector& lv2) {return lv1.E() > lv2.E(); };
  std::sort(BoostedDaughters->begin(), BoostedDaughters->end(), sortLambda);

  //  for(unsigned int i = 0; i < BoostedDaughters->size(); i++){
  //    std::cout << BoostedDaughters->at(i).E() << std::endl;
  //  }
  //  std::cout << std::endl;

  //Do the rotations
  float leadPhi = BoostedDaughters->begin()->Phi();
  float leadTheta = BoostedDaughters->begin()->Theta();
  float leadE =  BoostedDaughters->at(0).E();

  //  std::cout << leadPhi << ", " << BoostedDaughters->at(0).Phi() << std::endl;

  BoostedDaughters->at(0).RotateZ(-leadPhi);
  BoostedDaughters->at(0).SetPy(0);
  BoostedDaughters->at(0).RotateY(TMath::Pi()/2.0 - leadTheta);
  BoostedDaughters->at(0).SetPz(0);
  float subPsi = TMath::ATan2(BoostedDaughters->at(1).Py(), BoostedDaughters->at(1).Pz());

  for(auto lv = BoostedDaughters->begin(); lv != BoostedDaughters->end(); lv++){
    (*lv).RotateX(subPsi - TMath::Pi()/2.0);
  }
  BoostedDaughters->at(0).SetPz(0);
  BoostedDaughters->at(1).SetPz(0);

  //Reflect if necessary
  float rightSum = 0;
  float leftSum = 0;
  float topSum = 0;
  float botSum = 0;
  for(auto lv = BoostedDaughters->begin(); lv != BoostedDaughters->end(); lv++){
    if (lv->CosTheta() > 0){
      topSum+=lv->E();
    }
    else if (lv->CosTheta() < 0){
      botSum+=lv->E();
    }

    if (lv->Phi() > 0){
      rightSum+=lv->E();
    }
    else if (lv->Phi() < 0){
      leftSum+=lv->E();
    }
  }
  //Initialize image with 0's in all bins
  for (int nx = 0; nx < 31; nx++){
    for (int ny= 0; ny < 31; ny++){
      Image[nx][ny] = 0;
    }
  }
  //find the x and y coordinates in phi, theta binned space
  //Then fill Image with normalized energy
  
  for(auto lv = BoostedDaughters->begin(); lv != BoostedDaughters->end(); lv++){
    int x_bin = -1;
    int y_bin = -1;
    if (topSum >= botSum){
      x_bin = static_cast<int>(31*(lv->CosTheta() + 1)/(2.0));
      x_bin = x_bin%31;
    }
    else{
      x_bin = static_cast<int>(31*(-lv->CosTheta() + 1)/(2.0));
      x_bin = x_bin%31;
    }
    if (rightSum >= leftSum){
      y_bin = static_cast<int>(31*(lv->Phi() + TMath::Pi())/(2.0 * TMath::Pi()));
      y_bin = y_bin%31;
    }
    else{
      y_bin = static_cast<int>(31*(-lv->Phi() + TMath::Pi())/(2.0 * TMath::Pi()));
      y_bin = y_bin%31;
    }
    Image[x_bin][y_bin] += lv->E()/leadE ;    
  }
  delete BoostedDaughters;
}

bool checkKinematicsOfJets(const std::vector<pat::Jet> &jets, int length){
  bool check = true;
  for (int i = 0; i < length; i++){
    check = check && (jets.at(i).pt() > 400) && (fabs(jets.at(i).eta()) < 2.4);
  }
  return check;
}

bool checkLengthOfSubJets(const std::vector<pat::Jet> &jets, int length){
  bool check = true;
  for (int i = 0; i < length; i++){
    check = check && (jets.at(i).subjets("SoftDropPuppi").size() >= 2) && (jets.at(i).numberOfDaughters() > 2);
  }
  return check;
}
