#ifndef BEST_EVALUATION_H
#define BEST_EVALUATION_H

#include <Math/VectorUtil.h>
//#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
//#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
//#include "DataFormats/Common/interface/AssociationMap.h"
#include "VLQAnalysis/BESTAnalyzer/interface/CacheHandler.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


class BESTEvaluation {
 public:
  BESTEvaluation(const CacheHandler* cache)
    : isConfigured_(false),
    cache_(cache)
    {}
  void configure(const edm::ParameterSet&);
  std::vector<float> getPrediction(const float HImage[31][31], const float TImage[31][31], const float WImage[31][31], const float ZImage[31][31], const std::vector<float> BESTInputs);
  ~BESTEvaluation() {}


 private:
  static float getNormalValue(float value, float mean, float sigma){
    float result = (value-mean)/sigma;
    return result;
  }

  bool isConfigured_;
  const CacheHandler* cache_;
  std::string name_;
  edm::FileInPath path_;
  edm::FileInPath meanPath_;
  std::vector<float> means_;
  std::vector<float> sigmas_;
  std::vector<tensorflow::TensorShape> inputShapes_;
  std::vector<std::string> inputNames_;
  tensorflow::NamedTensorList inputTensors_;
  size_t kHiggs_;
  size_t kTop_;
  size_t kW_;
  size_t kZ_;
  size_t kBEST_;
  std::string outputName_;
  int NumBESTInputs_ = 94;

};
#endif
