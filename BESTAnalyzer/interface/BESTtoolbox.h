//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// BESTtoolbox.h --------------------------------------------------------------
//=================================================================================
// Header file containing functions for use with CMS EDAnalyzer and EDProducer ----
///////////////////////////////////////////////////////////////////////////////////

// make sure the functions are not declared more than once
#ifndef BESTtoolbox_H
#define BESTtoolbox_H

// include files
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "PhysicsTools/CandUtils/interface/Thrust.h"
#include "TMath.h"
#include "TLorentzVector.h"

// Fast Jet Include files
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include "fastjet/tools/Filter.hh"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ActiveAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>

///////////////////////////////////////////////////////////////////////////////////
// Functions ----------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////////

// calculate Legendre Polynomials
float LegendreP(float x, int order);

// calculate Fox Wolfram moments
int FWMoments(std::vector<TLorentzVector> particles, double (&outputs)[5] );

// get jet's constituents
void getJetDaughters(std::vector<reco::Candidate * > &daughtersOfJet, const pat::Jet &jet);

// store the jet variables
void storeJetVariables(std::map<std::string, float> &BESTVars, const pat::Jet &jet, const std::vector<reco::VertexCompositePtrCandidate> &secVertices);

// store the secondary vertex variables
void storeSecVertexVariables(std::map<std::string, float> &treeVars, std::map< std::string, std::vector<float> > &jetVecVars,
                             TLorentzVector jet, std::vector<reco::VertexCompositePtrCandidate> secVertices);

// store the rest frame variables
void storeRestFrameVariables(std::map<std::string, float> &BESTVars, std::vector<reco::Candidate *> daughtersOfJet,
                             const pat::Jet &jet, std::string frame, float mass);

// make rest frame z axis the boost axis
void pboost( TVector3 pbeam, TVector3 plab, TLorentzVector &pboo );

void initBESTVars(std::map<std::string, float> &BESMap, std::vector<std::string> listOfBESTVars);

std::vector<float> orderBESTVars(std::map<std::string, float> BESMap, std::vector<std::string> listOfBESTVars);

void prepareBoostedImage(const pat::Jet &jet, std::vector<reco::Candidate *> daughtersOfJet, float Image[31][31], const float mass);

#endif
