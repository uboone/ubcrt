///////////////////////////////////////////////////////////////////////
// Class:       CRTFileManager
// Plugin Type: service (art v2_05_01)
// File:        CRTFileManager_service.cc
//
// Generated at Thu Jan  3 14:09:49 2019 by Herbert Greenlee using cetskelgen
// from cetlib version v1_21_00.
//
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "CRTFileManager.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "boost/date_time/c_local_time_adjustor.hpp"
#include "IFDH_service.h"


namespace {

  // Local function to compare CRT hits according to time.

  bool CRTHitComp(const crt::CRTHit& crthit1, const crt::CRTHit& crthit2)
  {
    bool result = (crthit1.ts0_s < crthit2.ts0_s ||
		   (crthit1.ts0_s == crthit2.ts0_s && crthit1.ts0_ns < crthit2.ts0_ns));
    return result;
  }
}

// Constructor.

crt::CRTFileManager::CRTFileManager(fhicl::ParameterSet const & p, art::ActivityRegistry & areg) :
  fDebug(p.get<bool>("debug")),
  fMaxFiles(p.get<unsigned int>("maxFiles")),
  fCRTHitLabel(p.get<std::string>("CRTHitLabel")),
  fCRTVersion(p.get<std::string>("ubversion_CRTHits")),
  fCRTVersionTop(p.get<std::string>("ubversion_CRTHits_top", fCRTVersion)),
  fCRTMetadataDir(p.get<std::string>("CRTMetadataDir", std::string())),
  fSchema(p.get<std::string>("Schema", std::string()))
{
  setenv("TZ", "CST6CDT", 1);  // Fermilab time zone.
  tzset();

  // Maybe choose schema based on environment variable CRT_SCHEMA.

  if(fSchema == std::string() && getenv("CRT_SCHEMA") != 0 && *getenv("CRT_SCHEMA") != 0) {
    fSchema = std::string(getenv("CRT_SCHEMA"));
  }

  // Maybe use default schema based on following criteria.
  // 1. Xrootd version.
  // 2. Availability of valid proxy.
  // 3. Availability of /pnfs mount.
  //
  // Generallly, schema "root" is most preferred, then "file," then "https."

  if(fSchema == std::string()) {
    std::cout << "SAM schema was not specified as a fcl parameter.\n"
              << "We will choose a default schema."
              << std::endl;

    // Check xrootd version.

    std::string xrootd_version("v0");
    const char* s = getenv("XROOTD_VERSION");
    if(s != 0 && *s != 0)
      xrootd_version = std::string(s);
    std::cout << "Xrootd version = " << xrootd_version << std::endl;
    if(xrootd_version >= "v5_5")
      fSchema = "root";
    else {

      // Xrootd version is too low for tokens.
      // We can still use schema "root" if we have a valid proxy.

      int rc = system("voms-proxy-info -exists");
      if(rc == 0) {
        std::cout << "We have a valid proxy." << std::endl;
        fSchema = "root";
      }
      else {

        // No proxy, check availability of /pnfs, which allows to use schema "file."

        struct stat info;
        int ok = stat("/pnfs", &info);
        if(ok == 0 && (info.st_mode & S_IFDIR) != 0) {
          std::cout << "Filesystem /pnfs is mounted." << std::endl;
          fSchema = "file";
        }
        else {

          // No xrootd, proxy, or pnfs.
          // Use last resort schema https.
          
          std::cout << "No /pnfs mount." << std::endl;
          fSchema = "https";
        }
      }
    }
    std::cout << "Using default schema = " << fSchema << std::endl;
  }

  // Message.

  std::cout << "CRTFileManager service configured.\n"
	    << "  Max files = " << fMaxFiles  << "\n"
	    << "  CRT Hit Label = " << fCRTHitLabel  << "\n"
	    << "  CRT Version = " << fCRTVersion << "\n"
	    << "  CRT Top Version = " << fCRTVersionTop << "\n"
	    << "  CRT metadata directory = " << fCRTMetadataDir << "\n"
            << "  SAM schema = " << fSchema << std::endl;
}


// Destructor.

crt::CRTFileManager::~CRTFileManager()
{
  // Delete local copies of fetched files, if any.

  for(auto const& crt_event : fCRTEvents) {
    std::string class_name(crt_event.second->getTFile()->ClassName());
    if(class_name == std::string("TFile") ) {
      std::string local_filename(crt_event.second->getTFile()->GetName());
      if(fSchema != std::string("file") && local_filename.substr(0, 5) != std::string("/pnfs")) {
        std::cout << "Deleting local file " << local_filename << std::endl;
        remove(local_filename.c_str());
      }
    }
  }


}


// Fetch CRT metadata from extracted metadata directory.

void crt::CRTFileManager::prefetch(const std::vector<boost::posix_time::ptime>& posix_times) const
{
  // Can't do anything if there is no extracted crt metadata.

  if(fCRTMetadataDir == std::string())
    return;

  std::cout << "Prefech CRT metadata." << std::endl;

  // Loop over posix times -> subdirectories.

  for(auto const& posix_time : posix_times) {
    std::string subdir = boost::posix_time::to_iso_string(posix_time).substr(0,8);

    // Loop over files in this subdirectory

    std::string dirpath = fCRTMetadataDir + "/" + subdir;
    DIR* dir = opendir(dirpath.c_str());
    struct dirent* entry;
    if(dir != 0) {
      while((entry = readdir(dir)) != 0) {
	std::string filename = entry->d_name;
	if(filename != "." && filename != "..") {

	  // Read this file and construct CRTFileInfo object.

	  CRTFileInfo fileinfo;
	  std::cout << "Reading " << filename << std::endl;
	  std::string filepath = dirpath + "/" + filename;
	  std::ifstream mdfile(filepath);
	  std::string line;
	  if(mdfile.is_open()) {

	    // Read CRT binary file name.

	    getline(mdfile, line);
	    //std::cout << "CRT binary file = " << line << std::endl;
	    fileinfo.fFileName = line;
	    if(fCRTFiles.count(fileinfo.fFileName) > 0)
	      std::cout << "Already in cache." << std::endl;
	    else {

	      // Not in cache.
	      // Read start time.

	      getline(mdfile, line);
	      //std::cout << "Start time = " << line << std::endl;
	      boost::posix_time::ptime start_time_ptime =
		boost::date_time::parse_delimited_time<boost::posix_time::ptime>(line, 'T');
	      std::string start_time2 = boost::posix_time::to_iso_extended_string(start_time_ptime);
	      if(start_time2 != line) {
		//std::cout << "Converted start time = " << start_time2 << std::endl;
		throw cet::exception("CRTFileManager") << "Problem converting start time.";
	      }
	      fileinfo.fStartTime = start_time_ptime;

	      // Read end time.

	      getline(mdfile, line);
	      //std::cout << "End time = " << line << std::endl;
	      boost::posix_time::ptime end_time_ptime =
		boost::date_time::parse_delimited_time<boost::posix_time::ptime>(line, 'T');
	      std::string end_time2 = boost::posix_time::to_iso_extended_string(end_time_ptime);
	      if(end_time2 != line) {
		//std::cout << "Converted end time = " << end_time2 << std::endl;
		throw cet::exception("CRTFileManager") << "Problem converting end time.";
	      }
	      fileinfo.fEndTime = end_time_ptime;

	      // Remainder of file contains information about swizzled files.
	      // Loop over remaining lines.

	      std::string ubproject, ubstage, ubversion, swizzled_file;
	      while(!mdfile.eof()) {
		getline(mdfile, line);

		// Parse line into four words.

		ubproject = std::string();
		ubstage = std::string();
		ubversion = std::string();
		swizzled_file = std::string();
		std::istringstream ss(line);
		ss >> ubproject;
		ss >> ubstage;
		ss >> ubversion;
		ss >> swizzled_file;

		// Do the same metadata selection as if we were querying from sam.

		if(swizzled_file != std::string() &&
		   ((ubversion == fCRTVersionTop &&
		     (ubstage == "crt_swizzle1a" ||
		      ubstage == "crt_swizzle1b" ||
		      ubstage == "crt_swizzle1c")) ||
		    (ubversion == fCRTVersion &&
		     (ubstage == "crt_swizzle2" ||
		      ubstage == "crt_swizzle3" || 
		      ubstage == "crt_swizzle4")))) {
		  //std::cout << "Adding this CRT swizzled file to cache." << std::endl;
		  fileinfo.fSwizzled.push_back(swizzled_file);

		} // End metadata selection.
	      } // End loop over swizzled files.
	    } // End if in cache.

	    // Done reading file.

	    mdfile.close();

	    // Add this file info to cache.

	    fCRTFiles[fileinfo.fFileName] = fileinfo;

	  } // End if is open.
	} // End if not "." or ".."
      } // End loop over directory contents.
      closedir(dir);
    } 
  }
  return;
}


// Use sam or extracted CRT metadata to find swizzled CRT files that match the event time stamp.

std::vector<std::string> crt::CRTFileManager::findMatchingCRTFiles(art::Timestamp event_time) const
{
  // Result vector.

  std::vector<std::string> crtrootfiles;

  // Get ifdh service.

  art::ServiceHandle<ifdh_ns::IFDH> ifdh;

  // Convert the art time stamp into a string format that is understandable to sam.

  unsigned long evt_time_sec = event_time.timeHigh();
  unsigned long evt_time_nsec = event_time.timeLow();
  
  unsigned long time_tpc = evt_time_sec*1000000000 + evt_time_nsec;   // Nanoseconds
  unsigned long time_tpc1= time_tpc/1000;                             // Microseconds.

  // Lower bound of earliest possible start time of matching CRT binary file.
  // One day earlier than even time.

  unsigned long time_tpc0 = time_tpc1 - 86400000000;  // Microseconds.
  
  const char* tz = getenv("TZ");
  std::string tzs(tz);
  if (fDebug)
    std::cout<<"time-zone "<<tzs<<std::endl;
  
  if(tz != 0 && *tz != 0) {
    std::string tzs(tz);
    if(tzs != std::string("CST6CDT")) {
      // Timezone is wrong, throw exception.
      throw cet::exception("CRTFileManager") << "Wrong timezone: " << tzs;
    }
  }
  else {
    // Timezone is not set.  Throw exception.
    throw cet::exception("CRTFileManager") << "Timezone not set.";
  }
  
  boost::posix_time::ptime const time_epoch(boost::gregorian::date(1970, 1, 1));
  boost::posix_time::ptime this_event_time = 
    time_epoch + boost::posix_time::microseconds(time_tpc1);
  typedef boost::date_time::c_local_adjustor<boost::posix_time::ptime> local_adj;
  boost::posix_time::ptime this_event_localtime = local_adj::utc_to_local(this_event_time);
  std::cout << "Local time of event: " 
	    << boost::posix_time::to_iso_extended_string(this_event_localtime)<<std::endl;
  std::cout << "UTC time of event:   "
	    << boost::posix_time::to_iso_extended_string(this_event_time)<<std::endl;

  boost::posix_time::ptime this_event_time0 = 
    time_epoch + boost::posix_time::microseconds(time_tpc0);
  boost::posix_time::ptime this_event_localtime0 = local_adj::utc_to_local(this_event_time0);

  // Check cached files.
  // Make two attempts.
  // Call prefetch after the first attempt if we don't find enough matching files on
  // the first attempt.

  boost::posix_time::time_duration one_minute(0, 1, 0, 0);

  for(int itry=0; itry<2; ++itry) {

    crtrootfiles.erase(crtrootfiles.begin(), crtrootfiles.end());

    // Call prefetch on second attempt.

    if(itry > 0) {
      std::vector<boost::posix_time::ptime> posix_times;
      posix_times.push_back(this_event_localtime);
      posix_times.push_back(this_event_localtime0);
      prefetch(posix_times);
    }

    for(const auto& crtfileinfo : fCRTFiles) {
      const std::string& crtfile = crtfileinfo.first;    // crt swizzled file.
      const CRTFileInfo& fileinfo = crtfileinfo.second;  // crt binary file.
      if(this_event_localtime >= fileinfo.fStartTime + one_minute && 
	 this_event_localtime <= fileinfo.fEndTime - one_minute) {
	std::cout << "Found cached file " << crtfile << std::endl;
	for(auto const& crt_swizzled : fileinfo.fSwizzled)
	  crtrootfiles.push_back(crt_swizzled);
      }
    }

    // Break out of loop after first attempt if number of matching files is sufficient,
    // or if there is no extracted crt metadata to read for second attempt.

    std::cout << crtrootfiles.size() << " matching CRT binary files found in cache." << std::endl;

    if(crtrootfiles.size() >= 6 || fCRTMetadataDir == std::string())
      break;
  }
  if(crtrootfiles.size() < 6) {

    // Didn't match enough cached files.

    crtrootfiles.erase(crtrootfiles.begin(), crtrootfiles.end());

    // Query CRT binery files.

    std::string stringTime = boost::posix_time::to_iso_extended_string(this_event_localtime);
    stringTime = "'"+stringTime+"'";
    std::ostringstream dim;
    dim << "file_format " << "crt-binaryraw"
	<<" and file_type " << "data"
	<<" and start_time < " << stringTime 
	<< " and end_time > " << stringTime
	<< " and file_name " 
	<< "ProdRun" << boost::posix_time::to_iso_string(this_event_localtime0).substr(0,8) << "%,"
	<< "ProdRun" << boost::posix_time::to_iso_string(this_event_localtime).substr(0,8) << "%";
  
    if (fDebug)
      std::cout<<"dim = "<<dim.str()<<std::endl;
  
    // List those crtdaq files:
    std::vector< std::string > crtfiles = ifdh->translateConstraints(dim.str());
    if(fDebug)
      std::cout << "Found " << crtfiles.size() << " CRT binary files." << std::endl;

    // Get metadata of binary CRT files.

    for(auto const & crtfile : crtfiles) {
      if(fDebug)
	std::cout << "\n" << crtfile << std::endl;
      std::string md = ifdh->getMetadata(crtfile);

      // Extract CRT binary file start and end time from metadata as strings.

      std::string start_time;
      size_t n = md.find("Start Time: ");
      size_t m = 0;
      if(n < std::string::npos) {
	n += 12;
	m = md.find("+", n);
	if(m > n && m < std::string::npos) {
	  start_time = md.substr(n, m-n);
	  if(fDebug)
	    std::cout << "Start time = " << start_time << std::endl;
	}
      }
      std::string end_time;
      n = md.find("End Time: ");
      if(n < std::string::npos) {
	n += 10;
	m = md.find("+", n);
	if(m > n && m < std::string::npos) {
	  end_time = md.substr(n, m-n);
	  if(fDebug)
	    std::cout << "End time = " << end_time << std::endl;
	}
      }

      // Convert start and end times to boost ptime.

      boost::posix_time::ptime start_time_ptime =
	boost::date_time::parse_delimited_time<boost::posix_time::ptime>(start_time, 'T');
      boost::posix_time::ptime end_time_ptime =
	boost::date_time::parse_delimited_time<boost::posix_time::ptime>(end_time, 'T');

      // Make sure we closed the loop.

      std::string start_time2 = boost::posix_time::to_iso_extended_string(start_time_ptime);
      std::string end_time2 = boost::posix_time::to_iso_extended_string(end_time_ptime);
      if(fDebug) {
	std::cout << "Start ptime " << start_time2 << std::endl;
	std::cout << "End ptime " << end_time2 << std::endl;
      }
      if(start_time2 != start_time || end_time2 != end_time) {
	std::cout << "Start ptime " << start_time2 << std::endl;
	std::cout << "End ptime " << end_time2 << std::endl;
	throw cet::exception("CRTFileManager") << "Problem converting start and end time.";
      }
  
      // Query CRT swizzled files that are children of this CRT binary file.

      std::ostringstream dim1;
      dim1 << "file_format " << "artroot"
	   <<" and ((ub_project.version " << fCRTVersionTop
	   << " and ub_project.stage crt_swizzle1a,crt_swizzle1b,crt_swizzle1c)"
	   << " or (ub_project.version " << fCRTVersion
	   << " and ub_project.stage crt_swizzle2,crt_swizzle3,crt_swizzle4))"
	   << " and ischildof: (file_name " << crtfile
	   << " with availability physical )";

      if (fDebug)
	std::cout << "dim1 = " << dim1.str() << std::endl;
    
      std::vector< std::string > tmprootfiles = ifdh->translateConstraints(dim1.str());
    
      std::cout << "Found " << tmprootfiles.size()
		<< " daughters of " << crtfile << std::endl;
      if (tmprootfiles.size()>0) {
	for (const auto& artrootchild : tmprootfiles) {
	  crtrootfiles.push_back(artrootchild);
	  if(fDebug)
	    std::cout << artrootchild << std::endl;
	}
      }

      // Construct a CRTFileInfo struct and stash it in our cache, if not already there.

      if(fCRTFiles.count(crtfile) == 0) {
	std::cout << "Adding " << crtfile << " to file cache." << std::endl;
	CRTFileInfo fileinfo;
	fileinfo.fFileName = crtfile;
	fileinfo.fStartTime = start_time_ptime;
	fileinfo.fEndTime = end_time_ptime;
	fileinfo.fSwizzled = tmprootfiles;
	fCRTFiles[crtfile] = fileinfo;
	std::cout << "File metadata cache contains " << fCRTFiles.size() << " files." << std::endl;
      }
    }
  }
  std::sort(crtrootfiles.begin(), crtrootfiles.end());
  std::cout<<"\nNumber of matching CRT swizzled files: "<<crtrootfiles.size()<<std::endl;
  for(const auto& crt_swizzled : crtrootfiles)
    std::cout << crt_swizzled << std::endl;
  if (!crtrootfiles.size())
    std::cout << "\n\t CRTFileManager_module: No child CRT files found" << std::endl;

  // Throw exception if there are fewer than six CRT files.

  if(crtrootfiles.size() < 6) {
    throw cet::exception("CRTFileManager") << "Too few matching CRT files: " 
					   << crtrootfiles.size() << "\n";
  }

  // Done.

  return crtrootfiles;
}

// Return an open gallery file for the specified file.
// If we already have an open gallery file in the file cache, return that.
// Otherwise open a new file and add it to the cache.

gallery::Event& crt::CRTFileManager::openFile(std::string file_name)
{
  // Do we have an open gallery event for this file?
  // If yes, return it.

  for(auto& crtfile : fCRTEvents) {
    if(crtfile.first == file_name) {

      // Yes, reuse the open gallery event.

      std::cout << "Reusing open gallery event for file " << file_name << std::endl;
      return *(crtfile.second);
    }
  }

  // If we get here, we don't have a cached open file, so we will need to open one.

  // Make sure that the open file cache won't exceed the maximum size.

  if(fCRTEvents.size() >= fMaxFiles) {
    std::cout << "Closing file " << fCRTEvents.front().first << std::endl;

    // Delete local copy, if any.

    std::string class_name(fCRTEvents.front().second->getTFile()->ClassName());
    if(class_name == std::string("TFile") ) {
      std::string local_filename(fCRTEvents.front().second->getTFile()->GetName());
      if(fSchema != std::string("file") && local_filename.substr(0, 5) != std::string("/pnfs")) {
        std::cout << "Deleting local file " << local_filename << std::endl;
        remove(local_filename.c_str());
      }
    }
    fCRTEvents.pop_front();
  }

  // Get url of file.

  art::ServiceHandle<ifdh_ns::IFDH> ifdh;
  std::vector< std::string > urls;
  bool locate_ok = false;
  try {
    urls = ifdh->locateFile(file_name, fSchema);
    locate_ok = true;
  }
  catch(...) {
    locate_ok = false;
  }
  if(!locate_ok || urls.size() == 0) {
    std::cout << "No locations found for root file " << file_name << std::endl;
    throw cet::exception("CRTFileManager") << "Could not locate CRT file: " 
					   << file_name << "\n";
  }

  // Filter locations.
  // If there is a "production" location, use that.
  // Otherwise, use the first location.
  // After filtering, vector urls should contain one element with only the preferred url.

  int np = 0;
  if(urls.size() > 1) {
    for(size_t i=0; i<urls.size(); ++i) {
      if(urls[i].find("/production/") < std::string::npos) {
        np = i;
        break;
      }
    }
  }

  if(np > 0)
    urls[0] = urls[np];
  if(urls.size() > 1)
    urls.erase(urls.begin()+1, urls.end());
  std::cout<<"URL: " << urls[0] << std::endl;

  // Maybe copy url to local file.
  //
  // Calling fetchInput usually won't do anything if the url is an xrootd url.
  //
  // If ifdh decides it needs to copy the url, it uses the following destination
  // directories, in priority order.
  //
  // 1.  $IFDH_DATA_DIR
  // 2.  $_CONDOR_SCRATCH_DIR
  // 3.  $TMPDIR
  // 4.  /var/tmp (last resort).
  //
  // Don't ever let ifdh copy to /var/tmp, because that will likely crash the machine.
  //

  if(getenv("IFDH_DATA_DIR") == 0 &&
     getenv("_CONDOR_SCRATCH_DIR") == 0 &&
     getenv("TMPDIR") == 0) {
    std::cout << "Setting ifdh destination to current directory." << std::endl;
    putenv(const_cast<char*>("TMPDIR=."));
  }

  std::vector<std::string> local_filenames;

  // If the url is a file:// url, don't use ifdh fetch, but just strip off the initial "file://".

  if(urls[0].substr(0, 7) == std::string("file://"))
    local_filenames.push_back(urls[0].substr(7));
  else
    local_filenames.push_back(ifdh->fetchInput(urls[0]));
  std::cout << "Local filename = " << local_filenames[0] << std::endl;

  // gallery, when fed the list of the xrootd URLs, internally calls TFile::Open()
  // to open the file-list
  // In interactive mode, you have to get your proxy authenticated by issuing:
  // voms-proxy-init -noregen -rfc -voms fermilab:/fermilab/uboone/Role=Analysis
  // when you would like to launch a 'lar -c ... ... ...'
  // In batch mode, this step is automatically done
  //
  // Update 04/2025:
  //
  // Proxies are being disabled.  Ifdh may require a token for authentication.
  // A valid bearer token should always be available in batch jobs.
  // Tokens can be gotten interactively by this command:
  // htgettoken -a htvaultprod.fnal.gov -i <experiment>
    
  // Make several attempts to open file.

  int mtry = 5;
  int wait = 0;
  bool open_ok = false;
  std::unique_ptr<gallery::Event> crt_event;

  while(mtry-- > 0 && !open_ok) {
    if(wait != 0) {
      std::cout << "Waiting " << wait << " seconds." << std::endl;
      sleep(wait);
      wait *= 2;
    }
    else
      wait = 1;

    try {
      std::cout << "Open file." << std::endl;
      crt_event = std::unique_ptr<gallery::Event>(new gallery::Event(local_filenames));
      open_ok = true;
    }
    catch(...) {
      open_ok = false;
    }

    if(open_ok)
      std::cout << "Open succeeded." << std::endl;
    else {
      std::cout << "Open failed." << std::endl;
    }
  }

  if(open_ok) {
    std::cout<< "Opened the CRT root file from URL" << std::endl;
    fCRTEvents.emplace_back(file_name, std::move(crt_event));
    std::cout << "New [collection of] CRTEvent[s]. Its size is  "
	      << fCRTEvents.back().second->numberOfEventsInFile() << "." << std::endl;
  }
  else {
    std::cout << "Failed to open CRT root file URL." << std::endl;
    //throw cet::exception("CRTFileManager") << "Could not open CRT file: " 
    //				   << file_name << "\n";

    // Don't give art the chance to clean up.

    abort();
  }

  // Done.

  return *(fCRTEvents.back().second);
}

// Rewind gallery file to first event.
// Return true if success, false if fail.

bool crt::CRTFileManager::rewind(gallery::Event& event) const
{
  // Make several attempts to rewind/reopen file.

  int mtry = 5;
  int wait = 0;
  bool rewind_ok = false;

  while(mtry-- > 0 && !rewind_ok) {
    if(wait != 0) {
      std::cout << "Waiting " << wait << " seconds." << std::endl;
      sleep(wait);
      wait *= 2;
    }
    else
      wait = 1;

    try {
      std::cout << "Rewinding file." << std::endl;
      event.toBegin();
      rewind_ok = true;
    }
    catch(...) {
      rewind_ok = false;
    }

    if(rewind_ok)
      std::cout << "Rewind succeeded." << std::endl;
    else
      std::cout << "Rewind failed." << std::endl;
  }

  // Done.

  return rewind_ok;
}

// Get next nonempty filtered and sorted CRT hit collection starting from 
// current gallery event.
// Return empty collection if end of file is reached.

void crt::CRTFileManager::get_crt_hits(gallery::Event& event,
				       std::vector<crt::CRTHit>& output_hits) const
{
  // Make sure output hit container is initially empty.

  output_hits.erase(output_hits.begin(), output_hits.end());

  // Event loop.

  for(; !event.atEnd() && output_hits.size() == 0; event.next()) {

    // Get crt hit handle for the current event..

    gallery::ValidHandle<std::vector<crt::CRTHit> > h =
      event.getValidHandle< std::vector<crt::CRTHit> >(fCRTHitLabel);

    // Copy and filter hits into output hits collection.

    output_hits.reserve(h->size());
    for(const auto& crthit : *h) {

      // Filter out hits with bad times.

      if(crthit.ts0_s > 1300000000)
	output_hits.push_back(crthit);
    }
    if(output_hits.size() != 0)
      break;
  }

  // Sort hits into increaseing time order.

  std::sort(output_hits.begin(), output_hits.end(), CRTHitComp);
}


DEFINE_ART_SERVICE(crt::CRTFileManager)
