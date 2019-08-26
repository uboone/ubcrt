////////////////////////////////////////////////////////////////////////
// Class:       ExternalVtx
// File:        ExternalVtx_module.cc
//
// Generated at Wed June 26 2019 by Elena Gramellini elenag@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// declare the package with the CRT hit information.
#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"


// include track information
#include "lardataobj/RecoBase/Vertex.h"

// Declare the file for the associations to be made.
#include "lardata/Utilities/AssociationUtil.h"

// C++
#include <memory>
#include <array>
#include <fstream>
#include <string>
#include <sstream>
#include <string>

// ROOT
#include <TTree.h>

class ExternalVtx;

class ExternalVtx : public art::EDProducer
{
public:
  explicit ExternalVtx(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ExternalVtx(ExternalVtx const &) = delete;
  ExternalVtx(ExternalVtx &&) = delete;
  ExternalVtx &operator=(ExternalVtx const &) = delete;
  ExternalVtx &operator=(ExternalVtx &&) = delete;

  // Use the 'beginJob()' function to save information in the TTree.
  void beginJob() override;

  // Required functions.
  void produce(art::Event &e) override;
  void readExtVertices( std::vector<art::Ptr<recob::Vertex> > & vtxlist );

private:
  std::ifstream* fInputFile;
  std::string    fInputFileName; ///< Name of text file containing events to simulate
};




ExternalVtx::ExternalVtx(fhicl::ParameterSet const &p)
  : EDProducer{p}
  , fInputFile(0)
  , fInputFileName{p.get<std::string>("InputFileName")} 
// Initialize member data here.
{
  produces< std::vector<recob::Vertex>   >(); 
}

void ExternalVtx::beginJob()
{

}




void ExternalVtx::produce(art::Event &e)
{
  std::unique_ptr<std::vector<recob::Vertex> > VtxCol(new std::vector<recob::Vertex>);
  
  // check that the file is still good
  fInputFile = new std::ifstream(fInputFileName.c_str());
  if( !fInputFile->good() )
    throw cet::exception("External File Module") << "input text file "
                                                 << fInputFileName
                                                 << " cannot be read in produce().\n";


  double pos[3] = {};
  unsigned int this_run    = 0; 
  unsigned int this_subrun = 0; 
  unsigned int this_event  = 0; 
  int nVtx        = 0;

  // read in line to get event number and number of particles
  std::string oneLine;
  while ( std::getline(*fInputFile, oneLine))
    {
      std::istringstream inputLine;
      inputLine.str(oneLine);
      if (inputLine >> this_run >> this_subrun >> this_event >> nVtx) {
	for(unsigned short i = 0; i < nVtx; ++i) {
	  std::getline(*fInputFile, oneLine);
	  inputLine.clear();
	  inputLine.str(oneLine);
	  inputLine >> pos[0] >> pos[1] >> pos[2];
	  if (this_run    ==    e.run() ) {
	    if (this_subrun == e.subRun() ) {
	      if (this_event  ==  e.event() ) {
		recob::Vertex thisVtx(pos);
		VtxCol->push_back(thisVtx);
	      }// Check for Evt
	    }// Check for subrun
	  }// Check for run
	}// End read vtx lines
      }// End read event line
    }// End read file
  
  // Add the data products to the event.
  e.put(std::move(VtxCol));
} // End of the filter module


DEFINE_ART_MODULE(ExternalVtx)

//  LocalWords:  XYZ
