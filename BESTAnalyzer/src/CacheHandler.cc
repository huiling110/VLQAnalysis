/* 
Class to create an instance of a tensorflow session that loads the trained BEST network from a .pb file as a graph

Code mostly lifted from Devin Taylor's DeepDiTau project: https://github.com/dntaylor/DevTools-Ntuplizer/blob/106X/src/DeepCache.cc
*/

#include "VLQAnalysis/BESTAnalyzer/interface/CacheHandler.h"

//Handle all the loading in the constructor
CacheHandler::CacheHandler(const edm::FileInPath GraphPath){ //, const std::string GraphPath){
  tensorflow::SessionOptions options;
  tensorflow::setThreading(options, 1, "no_threads");
  std::string FullPath = GraphPath.fullPath();
  graph_.reset(tensorflow::loadGraphDef(FullPath)); //std::shared_ptr<tensorflow::GraphDef>
  session_ = tensorflow::createSession(graph_.get(), options);

}

CacheHandler::~CacheHandler(){
  tensorflow::closeSession(session_);
}
//Destructor jsut needs to close the TF session
