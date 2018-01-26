#include "uboone/CRT/CRTMerger.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "uboone/CRT/CRTProducts/MSetCRTFrag.hh"
#include "artdaq-core/Data/Fragments.hh"
#include "gallery/Event.h"
#include "CRTBernFEBDAQCore/Overlays/BernZMQFragment.hh"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "art/Framework/Principal/Handle.h"
#include "IFDH_service.h"
#include "uboone/RawData/utils/DAQHeaderTimeUBooNE.h"
#include "uboone/CRT/CRTProducts/CRTHit.hh"

#include <memory>
#include <string>
#include <vector>
#include <exception>
#include <sstream>
#include <unistd.h>
#include <ctime>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/local_time_adjustor.hpp"
#include <boost/date_time/c_local_time_adjustor.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
using namespace boost::posix_time;
crt::CRTMerger::CRTMerger(const fhicl::ParameterSet& pset): data_label_DAQHeader_(pset.get<std::string>("data_label_DAQHeader_"))
{
	std::cout<<"crt::CRTMerger::CRTMerger"<<std::endl;
	
	setenv("TZ", "CST+6CDT", 1);
	tzset();
	
	//setenv("IFDH_DATA_DIR","/uboone/data/users/kolahalb/MicroBooNE/",1);
	//std::cout<<"ifdh_data_dir "<<getenv("IFDH_DATA_DIR")<<std::endl;
	
	this->reconfigure(pset);
	produces< std::vector<crt::CRTHit> >();
	fTag		= pset.get<art::InputTag> ("InputTagName");
	_debug		= pset.get<bool>		 ("debug");
	fTimeOffSet	= pset.get<std::vector< unsigned long > > ("test_t_offset");
	previouscrtrootfile = "";
	if ( ! tIFDH ) tIFDH = new ifdh_ns::ifdh;
}

crt::CRTMerger::~CRTMerger()
{
	//fIFDH->cleanup();
	std::cout<<"crt::CRTMerger::~CRTmerger"<<std::endl;
	
	if ( ! tIFDH )
	delete tIFDH;
}

void crt::CRTMerger::produce(art::Event& event)
{
	if (_debug)
	{
		std::cout <<"crt::CRTMerger::produce, NEW EVENT" << std::endl;
	}
	
	//For this event
	if(_debug)
	std::cout<<fMaxCount<<std::endl;
	
	//get DAQ Header
	art::Handle< raw::DAQHeaderTimeUBooNE > rawHandle_DAQHeader;
	event.getByLabel(data_label_DAQHeader_, rawHandle_DAQHeader);
	
	if(!rawHandle_DAQHeader.isValid())
	{
		std::cout << "Run " << event.run() << ", subrun " << event.subRun()<< ", event " << event.event() << " has zero"<< " DAQHeaderTimeUBooNE  " << " in with label " << data_label_DAQHeader_ << std::endl;
    		return;
	}
	
	raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
	art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
	art::Timestamp evtTimeNTP = my_DAQHeader.ntp_time();
	
	std::cout<<"evt_timeGPS_sec "<<evtTimeGPS.timeHigh()<<",  evt_timeGPS_nsec "<<evtTimeGPS.timeLow()<<",  evt_timeNTP_sec "<<evtTimeNTP.timeHigh()<<",  evt_timeNTP_nsec "<<evtTimeNTP.timeLow()<<std::endl;
	
	// First find the art event time stamp
	//art::Timestamp evtTime = event.time();
	art::Timestamp evtTime = evtTimeGPS;
	
	unsigned long evt_time_sec = evtTime.timeHigh();//+fTimeOffSet[2];	
	unsigned long evt_time_nsec = evtTime.timeLow();
	
	// Use the information for configuring SAM query about the coincident crt-binrary-raw files
	unsigned long time_tpc = evt_time_sec*1000000000 + evt_time_nsec;
	unsigned long time_tpc1= time_tpc/1000;
	
	const char* tz = getenv("TZ");
	std::string tzs(tz);
	if (_debug)
	std::cout<<"time-zone "<<tzs<<std::endl;
	
	if(tz != 0 && *tz != 0)
	{
		std::string tzs(tz);
		if(tzs != std::string("CST+6CDT"))
		{
			// Timezone is wrong, throw exception.
			throw cet::exception("CRTMerger") << "Wrong timezone: " << tzs;
		}
	}
	else
	{
		// Timezone is not set.  Throw exception.
		throw cet::exception("CRTMerger") << "Timezone not set.";
	}
	
	boost::posix_time::ptime const time_epoch(boost::gregorian::date(1970, 1, 1));
	boost::posix_time::ptime this_event_time = time_epoch + boost::posix_time::microseconds(time_tpc1);
	typedef boost::date_time::c_local_adjustor<ptime> local_adj;
	ptime this_event_localtime = local_adj::utc_to_local(this_event_time);
	std::cout<<"Local time of event: "<< boost::posix_time::to_iso_extended_string(this_event_localtime)<<std::endl;
	
	if (_debug)
	std::cout << boost::posix_time::to_iso_extended_string(this_event_time)<<std::endl;
	
	std::string stringTime = boost::posix_time::to_iso_extended_string(this_event_localtime);
	stringTime = "'"+stringTime+"'";
	
	struct tm tm;
	memset(&tm, 0, sizeof(tm));
	tm.tm_isdst = -1;
	
	if (_debug)
	std::cout<<"TPC event time (s) "<<evt_time_sec<<"   event time (ns) "<<evt_time_nsec<<std::endl;
	
	// At this stage, find out the crtdaq files that may be time co-incident
	std::ostringstream dim;
	dim<<"file_format "<<"crt-binaryraw"<<" and file_type "<<"data"<<" and start_time < "<<stringTime<<" and end_time > "<<stringTime;
	
	if (_debug)
	std::cout<<"dim = "<<dim.str()<<std::endl;
	
	// List those crtdaq files:
	std::vector< std::string > crtfiles = tIFDH->translateConstraints(dim.str());
	
	// Find the corresponding child artroot file
	std::vector< std::string > crtrootfile;
	std::vector< std::string > tmprootfile;
	
	for(unsigned k =0; k<crtfiles.size(); k++)
	{
		std::ostringstream dim1;
		dim1<<"file_format "<<"artroot"<<" and ischildof: (file_name "<<crtfiles[k]<<")"<<std::endl;
		
		if (_debug)
		std::cout<<"dim1 = "<<dim1.str()<<std::endl;
		
		tmprootfile = tIFDH->translateConstraints(dim1.str());
		
		if (tmprootfile.size()>0)
		crtrootfile.push_back(tmprootfile[0]);
	}
	std::cout<<"total: "<<crtrootfile.size()<<std::endl;
	
	std::unique_ptr<std::vector<crt::CRTHit> > CRTHitEventsSet(new std::vector<crt::CRTHit>); //collection of CRTHits for this event
	
	for(unsigned crf_index = 0; crf_index < crtrootfile.size(); crf_index++)
	{
	if (_debug)
	std::cout<<"The child artroot file is "<<crtrootfile[crf_index]<<std::endl;
	
	// Read off the root file by streaming over internet. Use xrootd URL
	// This is alternative to use gsiftp URL. In that approach we use the URL to ifdh::fetchInput the crt file to local directory and keep it open as long as we make use of it
	//xrootd URL
	std::string schema = "root";
	std::vector< std::string > crtrootFile_xrootd_url = tIFDH->locateFile(crtrootfile[crf_index], schema);
	std::cout<<"xrootd URL: "<<crtrootFile_xrootd_url[0]<<std::endl;
	
	// gallery, when fed the list of the xrootd URL, internally calls TFile::Open() to open the file-list
	// In interactive mode, you have to get your proxy authenticated by issuing:
	// voms-proxy-init -noregen -rfc -voms fermilab:/fermilab/uboone/Role=Analysis
	// when you would like to launch a 'lar -c ... ... ...'
	// In batch mode, this step is automatically done
	
	gallery::Event fCRTEvent(crtrootFile_xrootd_url);
	std::cout<<"Opened the CRT root file from xrootd URL"<<std::endl;
	
	if (_debug)
	std::cout<<"Start merging attempts"<<std::endl;
	
	unsigned int merging = 0;
	//int triscuit = 0;
	for(fCRTEvent.toBegin(); !fCRTEvent.atEnd(); ++fCRTEvent)
	{
		//art::Handle< std::vector< crt::CRTHit> > rawHandle;
		//fCRTEvent.getByLabel(cTag, rawHandle);
		
		std::vector< crt::CRTHit >  const& CRTHitCollection = *(fCRTEvent.getValidHandle< std::vector<crt::CRTHit> >(cTag));	
		
		for(std::vector<int>::size_type hit_index=0; hit_index!=CRTHitCollection.size(); hit_index++)
		{
			crt::CRTHit CRTHitevent = CRTHitCollection[hit_index];
			unsigned long CRTtime_s	= CRTHitevent.ts0_s;
			unsigned long CRTtime_ns= CRTHitevent.ts0_ns;
			
			//unsigned long TPCtime_ns= evt_time_sec*1000000000+evt_time_nsec;
			unsigned long TPCtime_s	= evt_time_sec;
			unsigned long TPCtime_ns= evt_time_nsec;
			
			unsigned long MergingWindow_start = TPCtime_ns - 2000000;
			unsigned long MergingWindow_end	  = TPCtime_ns + 4000000;
			
			//std::cout<<"TPC_ns: "<<TPCtime_ns<<", CRT_ns: "<<CRTtime_ns<<", Diff "<<(TPCtime_ns-CRTtime_ns)/fTimeWindow<<std::endl;
			if(std::abs(CRTtime_s - TPCtime_s)<2)
			{
				if ((CRTtime_ns > MergingWindow_start) && (CRTtime_ns < MergingWindow_end))
				{
					std::cout<<"found match"<<std::endl;
				CRTHitEventsSet->emplace_back(CRTHitevent);
				merging += 1;
				//break;
			}
		}
	}
	//std::cout<<"# merging in the stream: "<<merging<<std::endl;
	/*
	if (merging==1)
	{
		std::cout<<"found match"<<std::endl;
		break;
	}
	*/
	//std::cout<<"Next CRT event ..."<<std::endl;
}
std::cout<<"# merging in the stream: "<<merging<<std::endl;
}
std::cout<<"# of CRTHitEventsSet before writing: "<<CRTHitEventsSet->size()<<std::endl;
	event.put(std::move(CRTHitEventsSet));
	
	if (_debug)
	std::cout<<"---X---"<<std::endl;
}

void crt::CRTMerger::reconfigure(fhicl::ParameterSet const & pset)
{
	std::cout<<"crt::CRTMerger::reconfigure"<<std::endl;
	
	cTag = {pset.get<std::string>("data_label_CRTHit_")};
	fTag = {pset.get<std::string>("InputTagName","crthit")};
	fTimeWindow = pset.get<unsigned>("TimeWindow",5000000);
}
DEFINE_ART_MODULE(crt::CRTMerger)
