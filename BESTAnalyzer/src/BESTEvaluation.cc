#include "VLQAnalysis/BESTAnalyzer/interface/BESTEvaluation.h"
#include <fstream>
#include <sstream>



void BESTEvaluation::configure(const edm::ParameterSet& iConfig){
  if(isConfigured_) return;
  name_ = iConfig.getParameter<std::string>("name");
  path_ = iConfig.getParameter<edm::FileInPath>("path");
  meanPath_ = iConfig.getParameter<edm::FileInPath>("means");

  std::string fullMeansPath = meanPath_.fullPath();
  std::ifstream inputMeansFile(fullMeansPath);
  std::string line, word;
  float mean, sigma;

  while (std::getline(inputMeansFile, line)) {
    std::istringstream ss(line);
    std::getline(ss,word,',');
    mean = std::stof(word);
    std::getline(ss,word,',');
    sigma = sqrt(std::stof(word));
    means_.push_back(mean);
    sigmas_.push_back(sigma);
  }
  inputMeansFile.close();


  inputNames_.push_back("input_1");
  inputShapes_.push_back(tensorflow::TensorShape{1, 31, 31, 1});
  kHiggs_ = 0; //Tracking the order of the image frames, must match the order the network was constructed with
  inputNames_.push_back("input_2");
  inputShapes_.push_back(tensorflow::TensorShape{1, 31, 31, 1});
  kTop_ = 1;
  inputNames_.push_back("input_3");
  inputShapes_.push_back(tensorflow::TensorShape{1, 31, 31, 1});
  kW_ = 2;
  inputNames_.push_back("input_4");
  inputShapes_.push_back(tensorflow::TensorShape{1, 31, 31, 1});
  kZ_ = 3;
  inputNames_.push_back("input_5");
  inputShapes_.push_back(tensorflow::TensorShape{1, NumBESTInputs_}); 
  kBEST_ = 4;
  outputName_ = "dense_20/Softmax";


  //Now set the internal tensor list to the size of your inputs

  inputTensors_.resize(inputShapes_.size());


  //Now make each element have the correct name and shape

  
  for (size_t i=0; i<inputShapes_.size(); i++) {
    inputTensors_[i] = tensorflow::NamedTensor(inputNames_[i], tensorflow::Tensor(tensorflow::DT_FLOAT, inputShapes_.at(i)));
  }
  //This class must be constructed after the global cache was created and passed to it
  //Get the graph from the cache

  auto graph = cache_->getGraph();


  //The rest of this is a sanity check
  for (size_t i=0; i<inputShapes_.size(); i++){
    const auto& name = graph.node(i).name(); 
    auto it = std::find(inputNames_.begin(), inputNames_.end(), name);
    if (it==inputNames_.end()) { //Check if input layer name is in the graph
      throw cms::Exception("BESTEvaluation")
	<< "Processing graph " << name_ << ".\n"
	<< "Unknown input name " << name;
    }

    const auto& shape = graph.node(i).attr().at("shape").shape();
    int j = std::distance(inputNames_.begin(),it); //Do inputs in correct order, even if they weren't declared in correct order
    for (int d=1; d<inputShapes_.at(j).dims(); d++) {
      if (shape.dim(d).size() != inputShapes_.at(j).dim_size(d)) {
	throw cms::Exception("BESTEvaluation")
	  << "Number of inputs in graph does not match those expected for " << name_ << ".\n"
	  << "Expected input " << j << " dim " << d << " = " << inputShapes_.at(j).dim_size(d) << "."
	  << " Found " << shape.dim(d).size() << ".";
      }
    }
  }


  const auto& outName = graph.node(graph.node_size() - 1).name();

  if (outName!=outputName_) {
    throw cms::Exception("BESTEvaluation")
      << "Processing graph " << name_ << ".\n"
      << "Unexpected output name. Expected " << outputName_ << " found " << outName << ".";
  }
  isConfigured_ = true;
}

std::vector<float> BESTEvaluation::getPrediction(const float HImage[31][31], const float TImage[31][31], const float WImage[31][31], const float ZImage[31][31], const std::vector<float> &BESTInputs){
  std::vector<tensorflow::Tensor> pred_vector; //vector of predictions to allow for evaluation multiple jets, but only using one at the moment
  tensorflow::Tensor prediction;
  std::vector<float> NNoutputs;
  
  inputTensors_.at(kHiggs_).second.flat<float>().setZero();
  inputTensors_.at(kTop_).second.flat<float>().setZero();
  inputTensors_.at(kW_).second.flat<float>().setZero();
  inputTensors_.at(kZ_).second.flat<float>().setZero();
  inputTensors_.at(kBEST_).second.flat<float>().setZero();

  //Setup the image input tensors
  for (int nx =0; nx < 31; nx++){
    for (int ny=0; ny < 31; ny++){
      (inputTensors_.at(kHiggs_).second).tensor<float, 4>()(0, nx, ny, 0) =  HImage[nx][ny];
      (inputTensors_.at(kTop_).second).tensor<float, 4>()(0, nx, ny, 0) =  TImage[nx][ny];
      (inputTensors_.at(kW_).second).tensor<float, 4>()(0, nx, ny, 0) =  WImage[nx][ny];
      (inputTensors_.at(kZ_).second).tensor<float, 4>()(0, nx, ny, 0) =  ZImage[nx][ny];
    }
  }
  //Setup the BES variable inputs
  for (int n=0; n < NumBESTInputs_; n++){
    inputTensors_.at(kBEST_).second.matrix<float>()(0, n) = getNormalValue(BESTInputs.at(n), means_[n], sigmas_[n]); //BESTInputs MUST be constructed in the correct order
  }

  //Evaluate
  tensorflow::run(&(cache_->getSession()),
                  inputTensors_,
                  {outputName_},
                  &pred_vector);


  prediction = tensorflow::Tensor(tensorflow::DT_FLOAT, {1, 6}); //6 here for the number of outputs
  for (int k = 0; k < 6; ++k) {
    const float pred = pred_vector[0].flat<float>()(k);
    if (!(pred >= 0 && pred <= 1)) {
      throw cms::Exception("BESTEvaluation")
	<< "invalid prediction = " << pred << " for pred_index = " << k;
    } 
    prediction.matrix<float>()(0, k) = pred;
  }
  for (int i = 0; i < 6; i++){
    NNoutputs.push_back(prediction.matrix<float>()(0,i)); //0,5 determined by number of outputs
  }
  return NNoutputs;
}
