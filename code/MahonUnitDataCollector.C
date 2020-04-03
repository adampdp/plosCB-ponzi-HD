#include "Lens.h"
#include "MahonUnitDataCollector.h"
#include "CG_MahonUnitDataCollector.h"
#include "NodeDescriptor.h"
#include "Node.h"
#include <map>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>



#define ITER getSimulation().getIteration()
#define RANK getSimulation().getRank()

void MahonUnitDataCollector::initialize(RNG& rng) 
{


  {std::ostringstream os1;
    os1<<"x1_"<<fileName;
    x1_file=new std::ofstream(os1.str().c_str(), std::ofstream::out | std::ofstream::app);}
  

  {std::ostringstream os1;
    os1<<"Spike_"<<fileName;
    spike_file=new std::ofstream(os1.str().c_str(), std::ofstream::out | std::ofstream::app);}


  {std::ostringstream os1;
    os1<<"x2_"<<fileName;
    x2_file=new std::ofstream(os1.str().c_str(), std::ofstream::out | std::ofstream::app);}
  

  {std::ostringstream os1;
    os1<<"x3_"<<fileName;
    x3_file=new std::ofstream(os1.str().c_str(), std::ofstream::out | std::ofstream::app);}


  assert(rows.size()==cols.size());
  assert(cols.size()==spikes.size());
  assert(cols.size()==x1.size());
  assert(cols.size()==x2.size());
  assert(cols.size()==x3.size());


}

void MahonUnitDataCollector::finalize(RNG& rng) 
{

  *spike_file<<std::endl;
  *x1_file<<std::endl;
  *x2_file<<std::endl;
  *x3_file<<std::endl;

  spike_file->close();
  x1_file->close();
  x2_file->close();
  x3_file->close();

  {std::ofstream ofs;
  ofs.open ("weights.dat", std::ofstream::out | std::ofstream::app);
  ofs << -10000 << " " << -10000 << " " << -10000 << std::endl;
  ofs.close();}


}

void MahonUnitDataCollector::dataCollectionLFP(Trigger* trigger, NDPairList* ndPairList) 
{

  std::ofstream& output=*x3_file;
  ShallowArray<double*>::iterator iter=x3.begin(), end=x3.end();
  double tot = 0;
  for (;iter!=end; ++iter) tot+=(**iter); 
   
  output<<ITER*deltaT<<" "<<tot<<std::endl;
  

}

void MahonUnitDataCollector::dataCollection(Trigger* trigger, NDPairList* ndPairList) 
{

  {std::ofstream& output=*x1_file;
  output<<ITER*deltaT<<" ";
  ShallowArray<double*>::iterator iter=x1.begin(), end=x1.end();
  for (int col=0; iter!=end && col<maxoutnum; ++iter, ++col) {
    output<<**iter<<" ";

  }
  output<<std::endl;}


 
  {std::ofstream& output=*x2_file;
  output<<ITER*deltaT<<" ";
  ShallowArray<double*>::iterator iter=x3.begin(), end=x3.end();
  for (int col=0; iter!=end && col<maxoutnum; ++iter, ++col) {
    output<<**iter<<" ";
   
  }
  output<<std::endl;}

}




void MahonUnitDataCollector::dataCollectionSpike(Trigger* trigger, NDPairList* ndPairList) 
{


  std::ofstream& output=*spike_file;
  
  ShallowArray<bool*>::iterator iter=spikes.begin(), end=spikes.end();
  for (int col=0; iter!=end; ++iter, ++col) {
    if (**iter == true) 
      output<<getSimulation().getIteration()<<" "<< col << std::endl;
    }

  

}




void MahonUnitDataCollector::getNodeIndices(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_MahonUnitDataCollectorInAttrPSet* CG_inAttrPset, CG_MahonUnitDataCollectorOutAttrPSet* CG_outAttrPset) 
{
  ShallowArray<unsigned,3,2> coords;
  CG_node->getNode()->getNodeCoords(coords);
  assert(coords.size()==2);
  rows.push_back(coords[0]);
  cols.push_back(coords[1]);
}

MahonUnitDataCollector::MahonUnitDataCollector() 
   : CG_MahonUnitDataCollector()
{
}

MahonUnitDataCollector::~MahonUnitDataCollector() 
{
  if (x1_file) delete x1_file; if (spike_file) delete spike_file;
  if (x2_file) delete x2_file; if (x3_file) delete x3_file; 
}

void MahonUnitDataCollector::duplicate(std::auto_ptr<MahonUnitDataCollector>& dup) const
{
   dup.reset(new MahonUnitDataCollector(*this));
}

void MahonUnitDataCollector::duplicate(std::auto_ptr<Variable>& dup) const
{
   dup.reset(new MahonUnitDataCollector(*this));
}

void MahonUnitDataCollector::duplicate(std::auto_ptr<CG_MahonUnitDataCollector>& dup) const
{
   dup.reset(new MahonUnitDataCollector(*this));
}

