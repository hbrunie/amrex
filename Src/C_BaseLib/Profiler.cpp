#ifdef BL_PROFILING

#include <Profiler.H>
#include <REAL.H>
#include <Utility.H>
#include <ParallelDescriptor.H>
#include <Array.H>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <stack>
#include <algorithm>
#include <limits>
#include <stdlib.h>


bool Profiler::bWriteAll = true;
bool Profiler::bNoOutput = false;
bool Profiler::bWriteFabs = true;
bool Profiler::bFirstCommWriteH = true;  // header
bool Profiler::bFirstCommWriteD = true;  // data
bool Profiler::bInitialized = false;

int Profiler::currentStep = 0;
int Profiler::csFlushSize = 2000000;
int Profiler::nProfFiles  = 64;
int Profiler::finestLevel = -1;
int Profiler::maxLevel    = -1;

Real Profiler::pctTimeLimit = 5.0;
Real Profiler::calcRunTime  = 0.0;
Real Profiler::startTime    = 0.0;
Real Profiler::timerTime    = 0.0;

Array<IntVect> Profiler::refRatio;
Array<Box> Profiler::probDomain;

std::stack<Real> Profiler::nestedTimeStack;
std::map<int, Real> Profiler::mStepMap;
std::map<std::string, Profiler::ProfStats> Profiler::mProfStats;
std::map<Real, std::string, std::greater<Real> > Profiler::mTimersTotalsSorted;
std::vector<Profiler::CommStats> Profiler::vCommStats;
std::map<std::string, Profiler *> Profiler::mFortProfs;
std::vector<std::string> Profiler::mFortProfsErrors;
const int mFortProfMaxErrors(32);
std::map<std::string, Profiler::CommFuncType> Profiler::CommStats::cftNames;
std::set<Profiler::CommFuncType> Profiler::CommStats::cftExclude;
int Profiler::CommStats::barrierNumber(0);
int Profiler::CommStats::reductionNumber(0);
int Profiler::CommStats::tagWrapNumber(0);
int Profiler::CommStats::tagMin(0);
int Profiler::CommStats::tagMax(0);
int Profiler::CommStats::csVersion(1);
std::vector<std::pair<std::string,int> > Profiler::CommStats::barrierNames;
std::vector<std::pair<int,int> > Profiler::CommStats::nameTags;
std::vector<std::string> Profiler::CommStats::nameTagNames;
std::vector<int> Profiler::CommStats::reductions;
std::vector<int> Profiler::CommStats::tagWraps;

std::string Profiler::procName("NoProcName");
int Profiler::procNumber(-1);
std::string Profiler::blProfDirName("bl_prof");
std::string Profiler::blCommProfDirName("bl_comm_prof");

#ifdef BL_CALL_TRACE
bool Profiler::bFirstTraceWriteH(true);  // header
bool Profiler::bFirstTraceWriteD(true);  // data
std::vector<Profiler::CallStats> Profiler::vCallTrace;
std::stack<int> Profiler::callIndexStack;
std::map<std::string, int> Profiler::mFNameNumbers;
int Profiler::callStackDepth(-1);
int Profiler::prevDepth(0);
std::map<std::string, int> Profiler::mRegionNameNumbers;
int Profiler::inNRegions(0);
std::vector<Profiler::RStartStop> Profiler::rStartStop;
int Profiler::CallStats::csVersion(1);
Real Profiler::CallStats::minCallTime(std::numeric_limits<Real>::max());
Real Profiler::CallStats::maxCallTime(std::numeric_limits<Real>::min());
const std::string Profiler::noRegionName("__NoRegion__");
#endif


Profiler::Profiler(const std::string &funcname)
    : bltstart(0.0), bltelapsed(0.0)
    , fname(funcname)
    , bRunning(false)
{
    start();
}


Profiler::~Profiler() {
  if(bRunning) {
    stop();
  }
}


void Profiler::Initialize() {
  if(bInitialized) {
    return;
  }
  int resultLen(-1);
  char cProcName[MPI_MAX_PROCESSOR_NAME + 11];
#ifdef BL_USE_MPI
  MPI_Get_processor_name(cProcName, &resultLen);
#endif
  if(resultLen < 1) {
    //strcpy(cProcName, "NoProcName");
    procName   = "NoProcName";
    procNumber = ParallelDescriptor::MyProc();
  } else {
    procName = cProcName;
#ifdef BL_HOPPER
    procNumber = atoi(procName.substr(3, std::string::npos).c_str());
#else
#ifdef BL_SIM_HOPPER
    //procNumber = (100 * ParallelDescriptor::MyProc()) % 6527;
    procNumber = ParallelDescriptor::MyProc();
    std::stringstream pname;
    pname << "nid" << procNumber;
    procName = pname.str();
    std::cout << ParallelDescriptor::MyProc() << "::procNumber = " << procNumber
              << "  procName = " << procName << std::endl;
#else
    procNumber = ParallelDescriptor::MyProc();
#endif
#endif
  }
  //std::cout << myProc << ":::: " << procName << "  len =  " << resultLen << std::endl;

  Real t0, t1;
  int nTimerTimes(1000);
  for(int i(0); i < nTimerTimes; ++i) {  // ---- time the timer
    t0 = ParallelDescriptor::second();
    t1 = ParallelDescriptor::second();
    timerTime += t1 - t0;
  }
  timerTime /= static_cast<Real> (nTimerTimes);

  startTime = ParallelDescriptor::second();

#ifdef BL_COMM_PROFILING
  vCommStats.reserve(csFlushSize);
#endif

#ifdef BL_CALL_TRACE
  vCallTrace.reserve(128000);
  BL_PROFILE_REGION_START(noRegionName);
#endif

  CommStats::cftExclude.insert(AllCFTypes);  // temporarily

  CommStats::cftNames["InvalidCFT"]     = InvalidCFT;
  CommStats::cftNames["AllReduceT"]     = AllReduceT; 
  CommStats::cftNames["AllReduceR"]     = AllReduceR; 
  CommStats::cftNames["AllReduceL"]     = AllReduceL; 
  CommStats::cftNames["AllReduceI"]     = AllReduceI; 
  CommStats::cftNames["AsendTsii"]      = AsendTsii;
  CommStats::cftNames["AsendTsiiM"]     = AsendTsiiM;
  CommStats::cftNames["AsendvTii"]      = AsendvTii;
  CommStats::cftNames["SendTsii"]       = SendTsii;
  CommStats::cftNames["SendvTii"]       = SendvTii;
  CommStats::cftNames["ArecvTsii"]      = ArecvTsii;
  CommStats::cftNames["ArecvTsiiM"]     = ArecvTsiiM;
  CommStats::cftNames["ArecvTii"]       = ArecvTii;
  CommStats::cftNames["ArecvvTii"]      = ArecvvTii;
  CommStats::cftNames["RecvTsii"]       = RecvTsii;
  CommStats::cftNames["RecvvTii"]       = RecvvTii;
  CommStats::cftNames["ReduceT"]        = ReduceT; 
  CommStats::cftNames["ReduceR"]        = ReduceR; 
  CommStats::cftNames["ReduceL"]        = ReduceL; 
  CommStats::cftNames["ReduceI"]        = ReduceI; 
  CommStats::cftNames["BCastTsi"]       = BCastTsi; 
  CommStats::cftNames["GatherTsT1Si"]   = GatherTsT1Si; 
  CommStats::cftNames["GatherTi"]       = GatherTi; 
  CommStats::cftNames["GatherRiRi"]     = GatherRiRi; 
  CommStats::cftNames["ScatterTsT1si"]  = ScatterTsT1si; 
  CommStats::cftNames["Barrier"]        = Barrier;
  CommStats::cftNames["Waitsome"]       = Waitsome;
  CommStats::cftNames["NameTag"]        = NameTag;
  CommStats::cftNames["AllCFTypes"]     = AllCFTypes;
  CommStats::cftNames["NoCFTypes"]      = NoCFTypes;
  CommStats::cftNames["IOStart"]        = IOStart;
  CommStats::cftNames["IOEnd"]          = IOEnd;
  CommStats::cftNames["TagWrap"]        = TagWrap;
  CommStats::cftNames["Allgather"]      = Allgather;
  CommStats::cftNames["Alltoall"]       = Alltoall;
  CommStats::cftNames["Alltoallv"]      = Alltoallv;
  CommStats::cftNames["Gatherv"]        = Gatherv;
  CommStats::cftNames["Get_count"]      = Get_count;
  CommStats::cftNames["Iprobe"]         = Iprobe;
  CommStats::cftNames["Test"]           = Test;
  CommStats::cftNames["Wait"]           = Wait;
  CommStats::cftNames["Waitall"]        = Waitall;

  // check for exclude file
  std::string exFile("CommFuncExclude.txt");
  std::vector<CommFuncType> vEx;

  Array<char> fileCharPtr;
  bool bExitOnError(false);  // in case the file does not exist
  ParallelDescriptor::ReadAndBcastFile(exFile, fileCharPtr, bExitOnError);

  CommStats::cftExclude.erase(AllCFTypes);

  if(fileCharPtr.size() > 0) {
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream cfex(fileCharPtrString, std::istringstream::in);

    while( ! cfex.eof()) {
        std::string cft;
        cfex >> cft;
        if( ! cfex.eof()) {
	  vEx.push_back(CommStats::StringToCFT(cft));
	}
    }
    for(int i(0); i < vEx.size(); ++i) {
      CommStats::cftExclude.insert(vEx[i]);
    }
  }
  bInitialized = true;
}


void Profiler::start() {
#ifdef _OPENMP
#pragma omp master
#endif
{
  ++mProfStats[fname].nCalls;
  bRunning = true;
  bltstart = ParallelDescriptor::second();
  nestedTimeStack.push(0.0);

#ifdef BL_CALL_TRACE
  int fnameNumber;
  std::map<std::string, int>::iterator it = Profiler::mFNameNumbers.find(fname);
  if(it == Profiler::mFNameNumbers.end()) {
    fnameNumber = Profiler::mFNameNumbers.size();
    Profiler::mFNameNumbers.insert(std::pair<std::string, int>(fname, fnameNumber));
  } else {
    fnameNumber = it->second;
  }
  ++callStackDepth;
  if(vCallTrace.size() == 0) {
    CallStats cs(callStackDepth, fnameNumber);
    ++cs.nCalls;
    cs.callTime = bltstart;
    vCallTrace.push_back(cs);
    CallStats::minCallTime = std::min(CallStats::minCallTime, bltstart);
    CallStats::maxCallTime = std::max(CallStats::maxCallTime, bltstart);
  } else {
    int topNameNumber(fnameNumber);
    if(vCallTrace.back().csFNameNumber == topNameNumber && callStackDepth != prevDepth) {
      ++(vCallTrace.back().nCalls);
    } else {
      CallStats cs(callStackDepth, fnameNumber);
      ++cs.nCalls;
      cs.callTime = bltstart;
      vCallTrace.push_back(cs);
      CallStats::minCallTime = std::min(CallStats::minCallTime, bltstart);
      CallStats::maxCallTime = std::max(CallStats::maxCallTime, bltstart);
    }
  }
  callIndexStack.push(vCallTrace.size() - 1);
  prevDepth = callStackDepth;
#endif
}
}

  
void Profiler::stop() {
#ifdef _OPENMP
#pragma omp master
#endif
{
  double tDiff(ParallelDescriptor::second() - bltstart);
  double nestedTime(0.0);
  bltelapsed += tDiff;
  bRunning = false;
  Real thisFuncTime(bltelapsed);
  if( ! nestedTimeStack.empty()) {
    nestedTime    = nestedTimeStack.top();
    thisFuncTime -= nestedTime;
    nestedTimeStack.pop();
  }
  if( ! nestedTimeStack.empty()) {
    nestedTimeStack.top() += bltelapsed;
  }
  mProfStats[fname].totalTime += thisFuncTime;

#ifdef BL_CALL_TRACE
  prevDepth = callStackDepth;
  --callStackDepth;
  if(vCallTrace.size() > 0) {
    if(vCallTrace.back().csFNameNumber == mFNameNumbers[fname]) {
      vCallTrace.back().totalTime = thisFuncTime + nestedTime;
      vCallTrace.back().stackTime = thisFuncTime;
    }
  }
  if( ! callIndexStack.empty()) {
    int index(callIndexStack.top());
    vCallTrace[index].totalTime = thisFuncTime + nestedTime;
    vCallTrace[index].stackTime = thisFuncTime;
    callIndexStack.pop();
  }
#endif
}
}


void Profiler::InitParams(const Real ptl, const bool writeall, const bool writefabs) {
  pctTimeLimit = ptl;
  bWriteAll = writeall;
  bWriteFabs = writefabs;
}


void Profiler::InitAMR(const int flev, const int mlev, const Array<IntVect> &rr,
                        const Array<Box> pd)
{
  finestLevel = flev;
  maxLevel    = mlev;
  refRatio.resize(rr.size());
  probDomain.resize(pd.size());
  for(int i(0); i < rr.size(); ++i) {
    refRatio[i] = rr[i];
  }
  for(int i(0); i < pd.size(); ++i) {
    probDomain[i] = pd[i];
  }
}


void Profiler::AddStep(const int snum) {
  currentStep = snum;
  mStepMap.insert(std::map<int, Real>::value_type(currentStep,
                                                  ParallelDescriptor::second()));
}


#ifdef BL_CALL_TRACE
void Profiler::RegionStart(const std::string &rname) {
  Real rsTime(ParallelDescriptor::second() - startTime);

  if(rname != noRegionName) {
    ++inNRegions;
  }
  if(inNRegions == 1) {
    RegionStop(noRegionName);
  }

  int rnameNumber;
  std::map<std::string, int>::iterator it = Profiler::mRegionNameNumbers.find(rname);
  if(it == Profiler::mRegionNameNumbers.end()) {
    rnameNumber = Profiler::mRegionNameNumbers.size();
    Profiler::mRegionNameNumbers.insert(std::pair<std::string, int>(rname, rnameNumber));
  } else {
    rnameNumber = it->second;
  }
  rStartStop.push_back(RStartStop(true, rnameNumber, rsTime));
}


void Profiler::RegionStop(const std::string &rname) {
  Real rsTime(ParallelDescriptor::second() - startTime);

  int rnameNumber;
  std::map<std::string, int>::iterator it = Profiler::mRegionNameNumbers.find(rname);
  if(it == Profiler::mRegionNameNumbers.end()) {  // ---- error
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "-------- error in RegionStop:  region " << rname
                << " never started."  << std::endl;
    }
    rnameNumber = Profiler::mRegionNameNumbers.size();
    Profiler::mRegionNameNumbers.insert(std::pair<std::string, int>(rname, rnameNumber));
  } else {
    rnameNumber = it->second;
  }
  rStartStop.push_back(RStartStop(false, rnameNumber, rsTime));

  if(rname != noRegionName) {
    --inNRegions;
  }
  if(inNRegions == 0) {
    RegionStart(noRegionName);
  }
}
#endif


void Profiler::Finalize() {
  if( ! bInitialized) {
    return;
  }
  if(bNoOutput) {
    bInitialized = false;
    return;
  }

  // --------------------------------------- gather global stats
  Real finalizeStart = ParallelDescriptor::second();  // time the timer
  const int nProcs(ParallelDescriptor::NProcs());
  const int myProc(ParallelDescriptor::MyProc());
  const int iopNum(ParallelDescriptor::IOProcessorNumber());

  // filter out profiler communications.
  CommStats::cftExclude.insert(AllCFTypes);

  int maxlen(0);
  Array<Real> gtimes(1);
  Array<long> ncalls(1);
  if(ParallelDescriptor::IOProcessor()) {
    gtimes.resize(nProcs);
    ncalls.resize(nProcs);
  }


  // -------- make sure the set of profiled functions is the same on all processors
  int pfStringsSize(0);
  std::ostringstream pfStrings;
  if(ParallelDescriptor::IOProcessor()) {
    for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
        it != mProfStats.end(); ++it)
    {
      pfStrings << it->first << '\n';
    }
    pfStrings << std::ends;
    pfStringsSize = pfStrings.str().size();
  }
  ParallelDescriptor::Bcast(&pfStringsSize, 1);

  char *pfChar = new char[pfStringsSize];
  if(ParallelDescriptor::IOProcessor()) {
    std::strcpy(pfChar, pfStrings.str().c_str());
  }
  ParallelDescriptor::Bcast(pfChar, pfStringsSize);

  if( ! ParallelDescriptor::IOProcessor()) {
    std::istringstream pfIn(pfChar);
    std::string pfName;
    while( ! pfIn.eof()) {
      pfIn >> pfName;
      if( ! pfIn.eof()) {
        std::map<std::string, ProfStats>::const_iterator it = mProfStats.find(pfName);
        if(it == mProfStats.end()) {
	  ProfStats ps;
          mProfStats.insert(std::pair<std::string, ProfStats>(pfName, ps));
	  //std::cout << myProc << ":  #### ProfName not found, inserting:  "
	            //<< pfName << std::endl;
        }
      }
    }
  }

  // ------- we really need to send names that are not on the ioproc
  // ------- to the ioproc but for now we will punt
  // ------- make a copy in case there are names not on the ioproc
  std::map<std::string, ProfStats> mProfStatsCopy;
  std::istringstream pfIn(pfChar);
  std::string pfName;
  while( ! pfIn.eof()) {
    pfIn >> pfName;
    if( ! pfIn.eof()) {
      std::map<std::string, ProfStats>::const_iterator it = mProfStats.find(pfName);
      if(it != mProfStats.end()) {
        mProfStatsCopy.insert(std::pair<std::string, ProfStats>(it->first, it->second));
      } else {
	std::cout << myProc << ":  #### Unknown ProfName:  =>" << pfName << "<=" << std::endl;
      }
    }
  }
  delete [] pfChar;

  // ------- now check for names not on the ioproc
  for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
      it != mProfStats.end(); ++it)
  {
    std::map<std::string, ProfStats>::const_iterator itcopy =
                                         mProfStatsCopy.find(it->first);
    //if(itcopy == mProfStatsCopy.end()) {
      //std::cout << myProc << ":  #### ProfName not on ioproc:  "
                //<< it->first << std::endl;
    //}
  }


  // ---------------------------------- now collect global data onto the ioproc
  for(std::map<std::string, ProfStats>::const_iterator it = mProfStatsCopy.begin();
      it != mProfStatsCopy.end(); ++it)
  {
    std::string profName(it->first);
    int pnLen(profName.size());
    maxlen = std::max(maxlen, pnLen);
    char cpn[profName.size() + 1];
    strcpy(cpn, profName.c_str());
    cpn[profName.size()] = '\0';

    ParallelDescriptor::Bcast(cpn, profName.size() + 1);
    ProfStats &pstats = mProfStatsCopy[profName];
    if(nProcs == 1) {
      gtimes[0] = pstats.totalTime;
      ncalls[0] = pstats.nCalls;
    } else {
      ParallelDescriptor::Gather(&pstats.totalTime, 1, gtimes.dataPtr(), 1, iopNum);
      ParallelDescriptor::Gather(&pstats.nCalls, 1, ncalls.dataPtr(), 1, iopNum);
    }
    Real tsum(0.0), tmin(gtimes[0]), tmax(gtimes[0]), tavg(0.0);
    long ncsum(0);
    if(ParallelDescriptor::IOProcessor()) {
      for(int i(0); i < gtimes.size(); ++i) {
        tsum += gtimes[i];
        tmin = std::min(tmin, gtimes[i]);
        tmax = std::max(tmax, gtimes[i]);
      }
      ProfStats &pstats = mProfStats[profName];  // not the copy on ioproc
      tavg = tsum / static_cast<Real> (gtimes.size());
      pstats.minTime = tmin;
      pstats.maxTime = tmax;
      pstats.avgTime = tavg;

      for(int i(0); i < ncalls.size(); ++i) {
        ncsum += ncalls[i];
      }
      // uncomment for reporting total calls summed over all procs
      //pstats.nCalls = ncsum;
    }
  }

  // --------------------------------------- print global stats to cout
  if(ParallelDescriptor::IOProcessor()) {
    bool bWriteAvg(true);
    if(nProcs == 1) {
      bWriteAvg = false;
    }
    WriteStats(std::cout, bWriteAvg);

#ifdef BL_CALL_TRACE
#ifdef BL_CALL_TRACE_HTML
    WriteHTML();
#endif
#endif
  }


  // --------------------------------------- print all procs stats to a file
  if(bWriteAll) {
    // --------------------- start nfiles block
    std::string cdir(blProfDirName);
    if(ParallelDescriptor::IOProcessor()) {
      if(BoxLib::FileExists(cdir)) {
        std::string newoldname(cdir + ".old." + BoxLib::UniqueString());
	std::cout << "Profiler::Finalize():  " << cdir
	          << " exists.  Renaming to:  " << newoldname << std::endl;
        std::rename(cdir.c_str(), newoldname.c_str());
      }
      if( ! BoxLib::UtilCreateDirectory(cdir, 0755)) {
        BoxLib::CreateDirectoryFailed(cdir);
      }
    }
    // Force other processors to wait until directory is built.
    ParallelDescriptor::Barrier("Profiler::Finalize::waitfordir");

    const int   myProc    = ParallelDescriptor::MyProc();
    const int   nProcs    = ParallelDescriptor::NProcs();
    const int   nOutFiles = std::max(1, std::min(nProcs, nProfFiles));
    const int   nSets     = (nProcs + (nOutFiles - 1)) / nOutFiles;
    const int   mySet     = myProc/nOutFiles;
    std::string cFileName(cdir + '/' + cdir + "_D_");
    std::string FullName  = BoxLib::Concatenate(cFileName, myProc % nOutFiles, 4);

    for(int iSet = 0; iSet < nSets; ++iSet) {
      if(mySet == iSet) {
        {  // scope
          std::ofstream csFile;

          if(iSet == 0) {   // First set.
            csFile.open(FullName.c_str(),
                        std::ios::out|std::ios::trunc|std::ios::binary);
          } else {
            csFile.open(FullName.c_str(),
                        std::ios::out|std::ios::app|std::ios::binary);
            csFile.seekp(0, std::ios::end);   // set to eof
          }
          if( ! csFile.good()) {
            BoxLib::FileOpenFailed(FullName);
          }

          // ----------------------------- write to file here
          WriteStats(csFile);
          // ----------------------------- end write to file here

          csFile.flush();
          csFile.close();
        }  // end scope

        int iBuff     = 0;
        int wakeUpPID = (myProc + nOutFiles);
        int tag       = (myProc % nOutFiles);
        if(wakeUpPID < nProcs) {
          ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
        }
      }
      if(mySet == (iSet + 1)) {   // Next set waits.
        int iBuff;
        int waitForPID = (myProc - nOutFiles);
        int tag        = (myProc % nOutFiles);
        ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
      }
    }
    // --------------------- end nfiles block

#ifdef BL_CALL_TRACE
    BL_PROFILE_REGION_STOP(noRegionName);
    WriteCallTrace();     // --------------------- write call trace data
#endif


    ParallelDescriptor::Barrier("Profiler::Finalize");
  }
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Profiler::Finalize():  time:  "   // time the timer
              << ParallelDescriptor::second() - finalizeStart << std::endl;
  }


#ifdef BL_COMM_PROFILING
  WriteCommStats();       // --------------------- write communication data
#endif

  // report any fortran errors.  should really check with all procs, just iop for now
  if(ParallelDescriptor::IOProcessor()) {
    if(Profiler::mFortProfs.size() > 0) {
      std::cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF FORTRAN PROFILING UNSTOPPED ERRORS" << std::endl;
      for(std::map<std::string, Profiler *>::iterator it = Profiler::mFortProfs.begin();
          it != Profiler::mFortProfs.end(); ++it)
      {
        std::cout << "FFFF function not stopped:  fname ptr = " << it->first
	          << "  ---->" << it->second << "<----" << std::endl;
      }
      std::cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF END FORTRAN PROFILING UNSTOPPED ERRORS" << std::endl;
    }
    if(Profiler::mFortProfsErrors.size() > 0) {
      std::cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF FORTRAN PROFILING ERRORS" << std::endl;
      if(Profiler::mFortProfsErrors.size() >= mFortProfMaxErrors) {
        std::cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF MAX FORTRAN ERRORS EXCEEDED" << std::endl;
      }
      for(int i(0); i < Profiler::mFortProfsErrors.size(); ++i) {
        std::cout << "FFFF " << Profiler::mFortProfsErrors[i] << std::endl;
      }
      std::cout << "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF END FORTRAN PROFILING ERRORS" << std::endl;
    }
  }

  bInitialized = false;
}


void Profiler::WriteStats(std::ostream &ios, bool bwriteavg) {
  const int myProc(ParallelDescriptor::MyProc());
  const int colWidth(10);

  mTimersTotalsSorted.clear();

  Real totalTimers(0.0), percent(0.0);
  int maxlen(0);
  for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
      it != mProfStats.end(); ++it)
  {
    std::string profName(it->first);
    int pnLen(profName.size());
    maxlen = std::max(maxlen, pnLen);

    if(bwriteavg) {
      totalTimers += it->second.avgTime;
    } else {
      totalTimers += it->second.totalTime;
    }
  }
  Real pTimeTotal(totalTimers);
  if(calcRunTime > 0.0 && bwriteavg == false) {
    pTimeTotal = calcRunTime;
  }

  ios << '\n' << '\n';
  if( ! bwriteavg) {
    ios << std::setfill('*')
        << std::setw(maxlen + 2 + 3 * (colWidth + 2) - (colWidth+12)) << "";
    ios << std::setfill(' ');
    ios << "  Processor:  " << std::setw(colWidth) << myProc << '\n';
  }

  // -------- write timers sorted by name
  WriteHeader(ios, colWidth, maxlen, bwriteavg);
  for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
      it != mProfStats.end(); ++it)
  {
    if(pTimeTotal > 0.0) {
      if(bwriteavg) {
        percent = 100.0 * (it->second.avgTime / pTimeTotal);
      } else {
        percent = 100.0 * (it->second.totalTime / pTimeTotal);
      }
    } else {
      percent = 100.0;
    }
    std::string fname(it->first);
    const ProfStats &pstats = it->second;
    WriteRow(ios, fname, pstats, percent, colWidth, maxlen, bwriteavg);
  }
  ios << '\n';
  ios << "Total Timers     = " << std::setw(colWidth) << totalTimers
      << " seconds." << '\n';
  if(calcRunTime > 0.0) {
    percent = 100.0 * totalTimers / calcRunTime;
    ios << "Calc Run Time    = " << std::setw(colWidth) << calcRunTime
        << " seconds." << '\n';
    ios << "Percent Coverage = " << std::setw(colWidth) << percent << " %" << '\n';
  }

  // -------- write timers sorted by percent
  ios << '\n' << '\n';
  WriteHeader(ios, colWidth, maxlen, bwriteavg);
  for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
      it != mProfStats.end(); ++it)
  {
    Real dsec;
    if(bwriteavg) {
      dsec = it->second.avgTime;
    } else {
      dsec = it->second.totalTime;
    }
    std::string sfir(it->first);
    mTimersTotalsSorted.insert(std::make_pair(dsec, sfir));
  }

  for(std::map<Real, std::string>::const_iterator it = mTimersTotalsSorted.begin();
      it != mTimersTotalsSorted.end(); ++it)
  {
    if(pTimeTotal > 0.0) {
      percent = 100.0 * (it->first / pTimeTotal);
    } else {
      percent = 100.0;
    }
    std::string fname(it->second);
    const ProfStats &pstats = mProfStats[fname];
    WriteRow(ios, fname, pstats, percent, colWidth, maxlen, bwriteavg);
  }
  if(bwriteavg) {
    ios << std::setfill('=') << std::setw(maxlen+4 + 5 * (colWidth+2)) << ""
        << '\n';
  } else {
    ios << std::setfill('=') << std::setw(maxlen+4 + 3 * (colWidth+2)) << ""
        << '\n';
  }
  ios << std::setfill(' ');
  ios << std::endl;


#ifdef BL_CALL_TRACE
  // -------- write timers sorted by inclusive times
  std::vector<std::string> fNumberNames(mFNameNumbers.size());
  for(std::map<std::string, int>::const_iterator it = Profiler::mFNameNumbers.begin();
      it != Profiler::mFNameNumbers.end(); ++it)
  {
    fNumberNames[it->second] = it->first;
  }
  ios << "vCallTrace.size() = " << vCallTrace.size() << std::endl;
  ios << "**************** ************vvvv" << '\n';
  for(int i(0); i < vCallTrace.size(); ++i) {
    CallStats &cs = vCallTrace[i];
    for(int indent(0); indent < cs.callStackDepth; ++indent) {
      ios << "----";
    }
    ios << "  " << fNumberNames[cs.csFNameNumber] << "  "
	<< cs.totalTime << "  " << cs.stackTime << '\n';
  }
  ios << "**************** ************^^^^" << '\n';

  // sort by total time
  std::vector<RIpair> funcTotalTimes(mFNameNumbers.size());
  for(int i(0); i < funcTotalTimes.size(); ++i) {
    funcTotalTimes[i].first  = 0.0;
    funcTotalTimes[i].second = i;
  }
  std::vector<int> callStack(1000, -1);  // use vector instead of stack for iterator
  for(int i(0); i < vCallTrace.size(); ++i) {
    CallStats &cs = vCallTrace[i];
    int depth(cs.callStackDepth);
    callStack[depth] = cs.csFNameNumber;
    bool recursiveCall(false);
    for(int d(0); d <  depth; ++d) {
      if(cs.csFNameNumber == callStack[d]) {
        ios << " RECURSIVE:  " << fNumberNames[cs.csFNameNumber] << '\n';
        recursiveCall = true;
      }
    }
    if( ! recursiveCall) {
      funcTotalTimes[cs.csFNameNumber].first += cs.totalTime;
    }
  }

  std::sort(funcTotalTimes.begin(), funcTotalTimes.end(), fTTComp());

  int numPrec(4);
  ios << '\n' << '\n';
  ios << std::setfill('-') << std::setw(maxlen+4 + 1 * (colWidth+2))
      << std::left << "Inclusive times " << '\n';
  ios << std::right << std::setfill(' ');
  ios << std::setw(maxlen + 2) << "Function Name"
      << std::setw(colWidth + 4) << "Time s"
      << '\n';

  for(int i(0); i < funcTotalTimes.size(); ++i) {
    ios << std::setw(maxlen + 2) << fNumberNames[funcTotalTimes[i].second] << "  "
        << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
        << funcTotalTimes[i].first << " s"
        << '\n';
  }
  ios << std::setfill('=') << std::setw(maxlen+4 + 1 * (colWidth+2)) << ""
      << '\n';
  ios << std::setfill(' ');
  ios << std::endl;

#endif
}


#ifdef BL_CALL_TRACE
void Profiler::WriteCallTrace(const bool bFlushing) {   // ---- write call trace data
    std::string cdir(blProfDirName);
    const int   myProc    = ParallelDescriptor::MyProc();
    const int   nProcs    = ParallelDescriptor::NProcs();
    const int   nOutFiles = std::max(1, std::min(nProcs, nProfFiles));
    const int   nSets     = (nProcs + (nOutFiles - 1)) / nOutFiles;
    const int   mySet     = myProc/nOutFiles;
    std::string cFilePrefix("bl_call_stats");
    std::string cFileName(cdir + '/' + cFilePrefix + "_D_");
    std::string FullName  = BoxLib::Concatenate(cFileName, myProc % nOutFiles, 4);

    if(ParallelDescriptor::IOProcessor()) {
      std::string globalHeaderFileName(cdir + '/' + cFilePrefix + "_H");
      std::ofstream csGlobalHeaderFile;
      csGlobalHeaderFile.open(globalHeaderFileName.c_str(),
                              std::ios::out | std::ios::trunc);
      if( ! csGlobalHeaderFile.good()) {
        BoxLib::FileOpenFailed(globalHeaderFileName);
      }
      csGlobalHeaderFile << "CallStatsProfVersion  " << CallStats::csVersion << '\n';
      csGlobalHeaderFile << "NProcs  " << nProcs << '\n';
      csGlobalHeaderFile << "NOutFiles  " << nOutFiles << '\n';

      for(std::map<std::string, int>::iterator it = mRegionNameNumbers.begin();
          it != mRegionNameNumbers.end(); ++it)
      {
        csGlobalHeaderFile << "RegionName " << '"' << it->first << '"'
	                   << ' ' << it->second << '\n';
      }
      for(int i(0); i < nOutFiles; ++i) {
        std::string headerName(cFilePrefix + "_H_");
        headerName = BoxLib::Concatenate(headerName, i, 4);
        csGlobalHeaderFile << "HeaderFile " << headerName << '\n';
      }
      csGlobalHeaderFile.flush();
      csGlobalHeaderFile.close();
    }

    std::string shortDFileName(cFilePrefix + "_D_");
    shortDFileName  = BoxLib::Concatenate(shortDFileName, myProc % nOutFiles, 4);
    std::string longDFileName(cdir + '/' + shortDFileName);

    std::string shortHeaderFileName(cFilePrefix + "_H_");
    shortHeaderFileName  = BoxLib::Concatenate(shortHeaderFileName, myProc % nOutFiles, 4);
    std::string longHeaderFileName(cdir + '/' + shortHeaderFileName);

    for(int iSet = 0; iSet < nSets; ++iSet) {
      if(mySet == iSet) {
        {  // scope
          std::ofstream csDFile, csHeaderFile;

          if(iSet == 0 && bFirstTraceWriteH) {   // First set.
	    bFirstTraceWriteH = false;
            csHeaderFile.open(longHeaderFileName.c_str(), std::ios::out | std::ios::trunc);
          } else {
            csHeaderFile.open(longHeaderFileName.c_str(), std::ios::out | std::ios::app);
            csHeaderFile.seekp(0, std::ios::end);   // set to eof
          }
          if( ! csHeaderFile.good()) {
            BoxLib::FileOpenFailed(longHeaderFileName);
          }

          if(iSet == 0 && bFirstTraceWriteD) {   // First set.
	    bFirstTraceWriteD = false;
            csDFile.open(longDFileName.c_str(),
	                 std::ios::out | std::ios::trunc | std::ios::binary);
          } else {
            csDFile.open(longDFileName.c_str(),
	                 std::ios::out | std::ios::app | std::ios::binary);
            csDFile.seekp(0, std::ios::end);   // set to eof
          }
          if( ! csDFile.good()) {
            BoxLib::FileOpenFailed(longDFileName);
          }

          // ----------------------------- write to file here
          csHeaderFile << "CallStatsProc " << myProc
		       << " nRSS " << rStartStop.size()
		       << " nTraceStats " << vCallTrace.size()
		       << "  datafile  " << shortDFileName
		       << "  seekpos  " << csDFile.tellp()
	               << '\n';
	  for(std::map<std::string, int>::iterator it = mFNameNumbers.begin();
	      it != mFNameNumbers.end(); ++it)
	  {
	    csHeaderFile << "fName " << '"' << it->first << '"'
	                 << ' ' << it->second << '\n';
	  }
	  csHeaderFile << std::setprecision(16) << "timeMinMax  "
	               << CallStats::minCallTime - startTime << ' '
	               << CallStats::maxCallTime - startTime << '\n';

          csHeaderFile.flush();
          csHeaderFile.close();

	  csDFile.write((char *) &rStartStop[0], rStartStop.size() * sizeof(RStartStop));
	  csDFile.write((char *) &vCallTrace[0], vCallTrace.size() * sizeof(CallStats));

	  csDFile.flush();
	  csDFile.close();
          // ----------------------------- end write to file here
        }  // end scope

        int iBuff     = 0;
        int wakeUpPID = (myProc + nOutFiles);
        int tag       = (myProc % nOutFiles);
        if(wakeUpPID < nProcs) {
          ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
        }
      }
      if(mySet == (iSet + 1)) {   // Next set waits.
        int iBuff;
        int waitForPID = (myProc - nOutFiles);
        int tag        = (myProc % nOutFiles);
        ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
      }
    }
    // --------------------- end nfiles block
}
#endif



#ifdef BL_CALL_TRACE
#ifdef BL_CALL_TRACE_HTML
void Profiler::WriteHTML() {
  std::vector<std::string> fNumberNames(mFNameNumbers.size());
  for(std::map<std::string, int>::const_iterator it = Profiler::mFNameNumbers.begin();
      it != Profiler::mFNameNumbers.end(); ++it)
  {
    fNumberNames[it->second] = it->first;
  }

  // write to html file
  std::string htmlFileName("CallTrace.html");
  std::stack<std::string> listEnds;
  std::ofstream csHTMLFile;
  csHTMLFile.open(htmlFileName.c_str(), std::ios::out | std::ios::trunc);
  if( ! csHTMLFile.good()) {
    BoxLib::FileOpenFailed(htmlFileName);
  }
  csHTMLFile << "<!DOCTYPE html>" << '\n';
  csHTMLFile << "<html>" << '\n';
  csHTMLFile << "<head>" << '\n';
  csHTMLFile << "<title>Call Tree</title>" << '\n';
  csHTMLFile << "</head>" << '\n';
  csHTMLFile << '\n';

  csHTMLFile << "<body>" << '\n';
  csHTMLFile << '\n';
  csHTMLFile << "<script type=\"text/javascript\">" << '\n';
  csHTMLFile << "function collapse(id) {" << '\n';
  csHTMLFile << "  var elem = document.getElementById(id);" << '\n';
  csHTMLFile << "  if(elem.style.display == '') {" << '\n';
  csHTMLFile << "    elem.style.display = 'none';" << '\n';
  csHTMLFile << "  } else {" << '\n';
  csHTMLFile << "    elem.style.display = '';" << '\n';
  csHTMLFile << "  }" << '\n';
  csHTMLFile << "}" << '\n';
  csHTMLFile << "</script>" << '\n';
  csHTMLFile << '\n';

  csHTMLFile << "<h3>Function calls.</h3>" << '\n';

  csHTMLFile << "<ul>" << '\n';
  listEnds.push("</ul>");

// the next two lines will indent the html
//#define IcsHTMLFile for(int id(0); id <= listEnds.size(); ++id) csHTMLFile << "  "; csHTMLFile
//#define IIcsHTMLFile for(int id(0); id < listEnds.size(); ++id) csHTMLFile << "  "; csHTMLFile
#define IcsHTMLFile csHTMLFile
#define IIcsHTMLFile csHTMLFile

  std::cout << "vCallTrace.size() = " << vCallTrace.size() << std::endl;
  for(int i(0); i < vCallTrace.size(); ++i) {
    CallStats &cs = vCallTrace[i];
    if(cs.nCalls > 1) {
      std::cout << "DDDDDDDDDD cs.nCalls = " << cs.nCalls << std::endl;
    }

    if(i == vCallTrace.size() - 1) {
        IcsHTMLFile << "<li>" << fNumberNames[cs.csFNameNumber] << "  "
	            << cs.totalTime << "  " << cs.stackTime
	            << "</li>" << '\n';
        for(int n(0); n < cs.callStackDepth; ++n) {
          IIcsHTMLFile << listEnds.top() << '\n';
          listEnds.pop();
          IIcsHTMLFile << listEnds.top() << '\n';
          listEnds.pop();
	}
    } else {
      CallStats &csNext = vCallTrace[i + 1];
      if(csNext.callStackDepth > cs.callStackDepth) {
        IcsHTMLFile << "<li>" << '\n';
        listEnds.push("</li>");
	IcsHTMLFile << "<a href=\"javascript:void(0)\" onclick=\"collapse('node" << i << "')\">"
	            << fNumberNames[cs.csFNameNumber] << "  "
	            << cs.totalTime << "  " << cs.stackTime
		    << "</a>" << '\n';
	if(cs.callStackDepth < 3) {
	  IcsHTMLFile << "<ul id=\"node" << i << "\" style=\"display:\">" << '\n';
	} else {
	  IcsHTMLFile << "<ul id=\"node" << i << "\" style=\"display:none\">" << '\n';
	}
        listEnds.push("</ul>");
      } else  if(csNext.callStackDepth == cs.callStackDepth) {
        IcsHTMLFile << "<li>" << fNumberNames[cs.csFNameNumber] << "  "
	            << cs.totalTime << "  " << cs.stackTime
	            << "</li>" << '\n';
      } else {
        IcsHTMLFile << "<li>" << fNumberNames[cs.csFNameNumber] << "  "
	            << cs.totalTime << "  " << cs.stackTime
	            << "</li>" << '\n';
        for(int n(0); n < cs.callStackDepth - csNext.callStackDepth; ++n) {
          IIcsHTMLFile << listEnds.top() << '\n';
          listEnds.pop();
          IIcsHTMLFile << listEnds.top() << '\n';
          listEnds.pop();
	}
      }
    }
  }

  if(listEnds.size() != 1) {
    std::cout << "**** Error:  listEnds.size() = " << listEnds.size() << std::endl;
  }
  csHTMLFile << listEnds.top() << '\n';
  listEnds.pop();

  csHTMLFile << "</body>" << '\n';
  csHTMLFile << "</html>" << '\n';

  csHTMLFile.close();
}
#endif
#endif


void Profiler::WriteCommStats(const bool bFlushing) {

  Real wcsStart(ParallelDescriptor::second());
  bool bAllCFTypesExcluded(OnExcludeList(AllCFTypes));
  if( ! bAllCFTypesExcluded) {
    CommStats::cftExclude.insert(AllCFTypes);  // temporarily
  }

  if(bFlushing) {
    int nCS(vCommStats.size());
    ParallelDescriptor::ReduceIntMax(nCS);
    if(nCS < csFlushSize) {
      if( ! bAllCFTypesExcluded) {
        CommStats::cftExclude.erase(AllCFTypes);
      }
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Bypassing flush, nCS < csFlushSize:  " << nCS
	          << "  " << csFlushSize << std::endl;
      }
      return;
    } else {
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Flushing commstats:  nCSmax = " << nCS << std::endl;
      }
    }
  }

  std::string cdir(blCommProfDirName);
  if(ParallelDescriptor::IOProcessor()) {
    if(bFirstCommWriteH) {
      if(BoxLib::FileExists(cdir)) {
        std::string newoldname(cdir + ".old." + BoxLib::UniqueString());
        std::cout << "Profiler::WriteCommStats():  " << cdir
                  << " exists.  Renaming to:  " << newoldname << std::endl;
        std::rename(cdir.c_str(), newoldname.c_str());
      }
    }
    if( ! BoxLib::UtilCreateDirectory(cdir, 0755)) {
      BoxLib::CreateDirectoryFailed(cdir);
    }
  }
  // Force other processors to wait until directory is built.
  ParallelDescriptor::Barrier("Profiler::WriteCommStats::waitfordir");

  bool bUseRelativeTimeStamp(true);
  if(bUseRelativeTimeStamp) {
    for(int ics(0); ics < Profiler::vCommStats.size(); ++ics) {
      CommStats &cs = Profiler::vCommStats[ics];
      cs.timeStamp -= startTime;
    }
  }


  // --------------------- start nfiles block
  const int   myProc    = ParallelDescriptor::MyProc();
  const int   nProcs    = ParallelDescriptor::NProcs();
  const int   nOutFiles = std::max(1, std::min(nProcs, nProfFiles));
  const int   nSets     = (nProcs + (nOutFiles - 1)) / nOutFiles;
  const int   mySet     = myProc/nOutFiles;

  std::string shortDFileName(cdir + "_D_");
  shortDFileName = BoxLib::Concatenate(shortDFileName, myProc % nOutFiles, 4);
  std::string longDFileName  = cdir + '/' + shortDFileName;

  std::string shortHeaderFileName(cdir + "_H_");
  shortHeaderFileName = BoxLib::Concatenate(shortHeaderFileName, myProc % nOutFiles, 4);
  std::string longHeaderFileName  = cdir + '/' + shortHeaderFileName;

  //double mpiWTick(MPI_Wtick());

  if(ParallelDescriptor::IOProcessor() && bFirstCommWriteH) {
    std::string globalHeaderFileName(cdir + '/' + cdir + "_H");
    std::ofstream csGlobalHeaderFile;
    csGlobalHeaderFile.open(globalHeaderFileName.c_str(), std::ios::out | std::ios::trunc);
    if( ! csGlobalHeaderFile.good()) {
      BoxLib::FileOpenFailed(globalHeaderFileName);
    }
    csGlobalHeaderFile << "CommProfVersion  " << CommStats::csVersion << '\n';
    csGlobalHeaderFile << "NProcs  " << nProcs << '\n';
    csGlobalHeaderFile << "CommStatsSize  " << sizeof(CommStats) << '\n';
    csGlobalHeaderFile << "NOutFiles  " << nOutFiles << '\n';
    csGlobalHeaderFile << "FinestLevel  " << finestLevel << '\n';
    csGlobalHeaderFile << "MaxLevel  " << maxLevel << '\n';
    for(int i(0); i < refRatio.size(); ++i) {
      csGlobalHeaderFile << "RefRatio  " << i << "  " << refRatio[i] << '\n';
    }
    for(int i(0); i < probDomain.size(); ++i) {
      csGlobalHeaderFile << "ProbDomain  " << i << "  " << probDomain[i] << '\n';
    }
    for(int i(0); i < nOutFiles; ++i) {
      std::string headerName(cdir + "_H_");
      headerName = BoxLib::Concatenate(headerName, i, 4);
      csGlobalHeaderFile << "HeaderFile " << headerName << '\n';
    }

    csGlobalHeaderFile.flush();
    csGlobalHeaderFile.close();
  }


  for(int iSet = 0; iSet < nSets; ++iSet) {
    if(mySet == iSet) {
      {  // scope
        std::ofstream csDFile, csHeaderFile;

        if(iSet == 0 && bFirstCommWriteH) {   // First set.
	  bFirstCommWriteH = false;
          csHeaderFile.open(longHeaderFileName.c_str(), std::ios::out | std::ios::trunc);
        } else {
          csHeaderFile.open(longHeaderFileName.c_str(), std::ios::out | std::ios::app);
          csHeaderFile.seekp(0, std::ios::end);   // set to eof
        }
        if( ! csHeaderFile.good()) {
          BoxLib::FileOpenFailed(longHeaderFileName);
        }

        if(iSet == 0 && bFirstCommWriteD) {   // First set.
	  bFirstCommWriteD = false;
          csDFile.open(longDFileName.c_str(),
                      std::ios::out | std::ios::trunc | std::ios::binary);
        } else {
          csDFile.open(longDFileName.c_str(),
                      std::ios::out | std::ios::app | std::ios::binary);
          csDFile.seekp(0, std::ios::end);   // set to eof
        }
        if( ! csDFile.good()) {
          BoxLib::FileOpenFailed(longDFileName);
        }

        // ----------------------------- write to file here
        csHeaderFile << "CommProfProc  " << myProc
                     << "  nCommStats  " << vCommStats.size()
                     << "  datafile  " << shortDFileName
	             << "  seekpos  " << csDFile.tellp()
		     << "  " << procName << '\n';
        for(int ib(0); ib < CommStats::barrierNames.size(); ++ib) {
          int seekindex(CommStats::barrierNames[ib].second);
          CommStats &cs = vCommStats[seekindex];
          csHeaderFile << "bNum  " << cs.tag  // tag is used for barrier number
                       << ' ' << '"' << CommStats::barrierNames[ib].first << '"'
                       << ' ' << seekindex << '\n';
        }
        for(int ib(0); ib < CommStats::nameTags.size(); ++ib) {
          int seekindex(CommStats::nameTags[ib].second);
          csHeaderFile << "nTag  " << CommStats::nameTags[ib].first << ' '
                       << seekindex << '\n';
        }
        for(int ib(0); ib < CommStats::reductions.size(); ++ib) {
          int seekindex(CommStats::reductions[ib]);
          CommStats &cs = vCommStats[seekindex];
          csHeaderFile << "red  " << cs.tag  // tag is used for reduction number
	               << ' ' << seekindex << '\n';
        }
	if(vCommStats.size() > 0) {
	  csHeaderFile << std::setprecision(16)
	               << "timeMinMax  " << vCommStats[0].timeStamp << ' '
	               << vCommStats[vCommStats.size()-1].timeStamp << '\n';
	  csHeaderFile << std::setprecision(16)
	               << "timerTime  " << timerTime << '\n';
	} else {
	  csHeaderFile << "timeMinMax  0.0  0.0" << '\n';
	}
        for(int i(0); i < CommStats::nameTagNames.size(); ++i) {
          csHeaderFile << "nameTagNames  " << '"' << CommStats::nameTagNames[i]
                       << '"' << '\n';
        }
        csHeaderFile << "tagRange  " << CommStats::tagMin << ' '
	             << CommStats::tagMax << '\n';
        for(int i(0); i < CommStats::tagWraps.size(); ++i) {
          csHeaderFile << "tagWraps  " << CommStats::tagWraps[i] << '\n';
        }
	csHeaderFile.flush();
        csHeaderFile.close();

	csDFile.write((char *) &vCommStats[0], vCommStats.size() * sizeof(CommStats));

        csDFile.flush();
        csDFile.close();
        // ----------------------------- end write to file here
      }  // end scope

      int iBuff     = 0;
      int wakeUpPID = (myProc + nOutFiles);
      int tag       = (myProc % nOutFiles);
      if(wakeUpPID < nProcs) {
        ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
      }
    }
    if(mySet == (iSet + 1)) {   // Next set waits.
      int iBuff;
      int waitForPID = (myProc - nOutFiles);
      int tag        = (myProc % nOutFiles);
      ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
    }
  }
  // --------------------- end nfiles block


  // --------------------- flush the data
  vCommStats.clear();
  CommStats::barrierNames.clear();
  CommStats::nameTags.clear();
  CommStats::reductions.clear();
  if( ! bAllCFTypesExcluded) {
    CommStats::cftExclude.erase(AllCFTypes);
  }

  ParallelDescriptor::Barrier("Profiler::WriteCommStats::end");

  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Profiler::WriteCommStats():  time:  "
              << ParallelDescriptor::second() - wcsStart << std::endl;
  }
}



void Profiler::WriteHeader(std::ostream &ios, const int colWidth,
                           const Real maxlen, const bool bwriteavg)
{
  if(bwriteavg) {
    ios << std::setfill('-') << std::setw(maxlen+4 + 5 * (colWidth+2))
        << std::left << "Total times " << '\n';
    ios << std::right << std::setfill(' ');
    ios << std::setw(maxlen + 2) << "Function Name"
        << std::setw(colWidth + 2) << "NCalls"
        << std::setw(colWidth + 2) << "Min"
        << std::setw(colWidth + 2) << "Avg"
        << std::setw(colWidth + 2) << "Max"
        << std::setw(colWidth + 4) << "Percent %"
        << '\n';
  } else {
    ios << std::setfill('-') << std::setw(maxlen+4 + 3 * (colWidth+2))
        << std::left << "Total times " << '\n';
    ios << std::right << std::setfill(' ');
    ios << std::setw(maxlen + 2) << "Function Name"
        << std::setw(colWidth + 2) << "NCalls"
        << std::setw(colWidth + 2) << "Time"
        << std::setw(colWidth + 4) << "Percent %"
        << '\n';
  }
}

void Profiler::WriteRow(std::ostream &ios, const std::string &fname,
                        const ProfStats &pstats, const Real percent,
			const int colWidth, const Real maxlen,
			const bool bwriteavg)
{
    int numPrec(4), pctPrec(2);
    if(bwriteavg) {
      ios << std::right;
      ios << std::setw(maxlen + 2) << fname << "  "
          << std::setw(colWidth) << pstats.nCalls << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << pstats.minTime << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << pstats.avgTime << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << pstats.maxTime << "  "
          << std::setprecision(pctPrec) << std::fixed << std::setw(colWidth)
	  << percent << " %" << '\n';
    } else {
      ios << std::setw(maxlen + 2) << fname << "  "
          << std::setw(colWidth) << pstats.nCalls << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << pstats.totalTime << "  "
          << std::setprecision(pctPrec) << std::fixed << std::setw(colWidth)
	  << percent << " %" << '\n';
    }
}


bool Profiler::OnExcludeList(CommFuncType cft) {
  // 
  // the idea for NoCFTypes is to allow local filtering/unfiltering
  // while preserving the users exclude list
  // possibly use a filter stack instead
  // might need caching if performance is a problem
  // what to do if both all and none are on the list?
  // 
  if(CommStats::cftExclude.empty()) {  // nothing on the exclude list
    return false;
  }
  std::set<CommFuncType>::iterator cfti;
  cfti = CommStats::cftExclude.find(NoCFTypes);
  if(cfti != CommStats::cftExclude.end()) {  // found, exclude nothing
    return false;
  }
  cfti = CommStats::cftExclude.find(cft);
  if(cfti != CommStats::cftExclude.end()) {  // found, exclude this cft
    return true;
  }
  cfti = CommStats::cftExclude.find(AllCFTypes);
  if(cfti != CommStats::cftExclude.end()) {  // found, exclude all types
    return true;
  }
  return false;
}


void Profiler::AddCommStat(const CommFuncType cft, const int size,
                           const int pid, const int tag)
{
  if(OnExcludeList(cft)) {
    return;
  }
  vCommStats.push_back(CommStats(cft, size, pid, tag, ParallelDescriptor::second()));
}


void Profiler::AddBarrier(const std::string &message, const bool beforecall) {
  const CommFuncType cft(Profiler::Barrier);
  if(OnExcludeList(cft)) {
    return;
  }
  if(beforecall) {
    int tag(CommStats::barrierNumber);
    vCommStats.push_back(CommStats(cft, 0, BeforeCall(), tag,
                                   ParallelDescriptor::second()));
    CommStats::barrierNames.push_back(std::make_pair(message, vCommStats.size() - 1));
    ++CommStats::barrierNumber;
  } else {
    int tag(CommStats::barrierNumber - 1);  // it was incremented before the call
    vCommStats.push_back(CommStats(cft, AfterCall(), AfterCall(), tag,
                                   ParallelDescriptor::second()));
  }
}


void Profiler::TagRange(const int mintag, const int maxtag) {
  CommStats::tagMin = mintag;
  CommStats::tagMax = maxtag;
}


void Profiler::AddTagWrap() {
  const CommFuncType cft(Profiler::TagWrap);
  if(OnExcludeList(cft)) {
    return;
  }
  ++CommStats::tagWrapNumber;
  int tag(CommStats::tagWrapNumber);
  int index(CommStats::nameTags.size());
  CommStats::tagWraps.push_back(index);
  vCommStats.push_back(CommStats(cft, index,  vCommStats.size(), tag,
                       ParallelDescriptor::second()));
}


void Profiler::AddAllReduce(const CommFuncType cft, const int size,
                            const bool beforecall)
{
  if(OnExcludeList(cft)) {
    return;
  }
  if(beforecall) {
    int tag(CommStats::reductionNumber);
    vCommStats.push_back(CommStats(cft, size, BeforeCall(), tag,
                                   ParallelDescriptor::second()));
    CommStats::reductions.push_back(vCommStats.size() - 1);
    ++CommStats::reductionNumber;
  } else {
    int tag(CommStats::reductionNumber - 1);
    vCommStats.push_back(CommStats(cft, size, AfterCall(), tag,
                                   ParallelDescriptor::second()));
  }
}


void Profiler::AddWaitsome(const CommFuncType cft, const Array<MPI_Request> &reqs,
                           const int completed, const Array<int> &indx,
			   const Array<MPI_Status> &status, const bool beforecall)
{
#ifdef BL_USE_MPI
  if(OnExcludeList(cft)) {
    return;
  }
  if(beforecall) {
    vCommStats.push_back(CommStats(cft, BeforeCall(), BeforeCall(), NoTag(),
                         ParallelDescriptor::second()));
  } else {
    for(int i(0); i < completed; ++i) {
      MPI_Status stat(status[i]);
      int c;
      BL_MPI_REQUIRE( MPI_Get_count(&stat, MPI_UNSIGNED_CHAR, &c) );
      vCommStats.push_back(CommStats(cft, c, stat.MPI_SOURCE, stat.MPI_TAG,
                           ParallelDescriptor::second()));
    }
  }
#endif
}


int Profiler::NameTagNameIndex(const std::string &name) {  // prob need to opt this
  for(int i(0); i < CommStats::nameTagNames.size(); ++i) {
    if(CommStats::nameTagNames[i] == name) {
      return i;
    }
  }
  CommStats::nameTagNames.push_back(name);
  return CommStats::nameTagNames.size() - 1;
}


void Profiler::AddNameTag(const std::string &name) {
  const CommFuncType cft(Profiler::NameTag);
  if(OnExcludeList(cft)) {
    return;
  }
  int tag(NameTagNameIndex(name));
  int index(CommStats::nameTags.size());
  vCommStats.push_back(CommStats(cft, index,  vCommStats.size(), tag,
                       ParallelDescriptor::second()));
  CommStats::nameTags.push_back(std::make_pair(tag, vCommStats.size() - 1));
}


Profiler::CommFuncType Profiler::CommStats::StringToCFT(const std::string &s) {
  return CommStats::cftNames[s];
}


void Profiler::CommStats::Filter(CommFuncType cft) {
  if( ! OnExcludeList(cft)) {
    CommStats::cftExclude.insert(cft);
  }
}


void Profiler::CommStats::UnFilter(CommFuncType cft) {
  if(OnExcludeList(cft)) {
    CommStats::cftExclude.erase(cft);
  }
}


void Profiler::PerfMonProcess() {
#ifdef BL_USE_MPI
  MPI_Status status;
  bool finished(false);
  int recstep(-1), rtag(0);
  while( ! finished) {
    std::cout << "**** _in PerfMonProcess:  waiting for rtag = " << rtag << std::endl;
    MPI_Recv(&recstep, 1, MPI_INT, 0, rtag, ParallelDescriptor::CommunicatorInter(), &status);
    std::cout << "**** _in PerfMonProcess:  recv step = " << recstep << std::endl;
    ++rtag;
    if(recstep < 0) {
      finished = true;
    }
  }
#endif
  std::cout << "**** _in PerfMonProcess:  exiting." << std::endl;
}



namespace {
  const int EOS(-1);

  std::string Trim(const std::string &str) {
    int n;
    for(n = str.size(); --n >= 0; ) {
      if(str[n] != ' ' ) {
	break;
      }
    }
    std::string result;
    for(int i(0); i <= n; ++i) {
      result += str[i];
    }
    return result;
  }

  std::string Fint_2_string(const int *iarr, int nlen) {
    std::string res;
    for(int i(0); i < nlen && *iarr != EOS; ++i) {
      res += *iarr++;
    }
    return Trim(res);
  }
}

#include <BLFort.H>

BL_FORT_PROC_DECL(BL_PROFFORTFUNCSTART_CPP, bl_proffortfuncstart_cpp)
  (const int istr[], const int *NSTR)
{
  std::string fName(Fint_2_string(istr, *NSTR));
  std::map<std::string, Profiler *>::const_iterator it = Profiler::mFortProfs.find(fName);
  if(it == Profiler::mFortProfs.end()) {  // make a new profiler
    Profiler::mFortProfs.insert(std::pair<std::string, Profiler *>(fName, new Profiler(fName)));
  } else {  // error:  fname is already being profiled
    std::string estring("bl_proffortfuncstart error:  mFortProfs function already being profiled:  ");
    estring += fName;
    if(Profiler::mFortProfsErrors.size() < mFortProfMaxErrors) {
      Profiler::mFortProfsErrors.push_back(estring);
    }
  }
}

BL_FORT_PROC_DECL(BL_PROFFORTFUNCSTOP_CPP, bl_proffortfuncstop_cpp)
  (const int istr[], const int *NSTR)
{
  std::string fName(Fint_2_string(istr, *NSTR));
  std::map<std::string, Profiler *>::const_iterator it = Profiler::mFortProfs.find(fName);
  if(it == Profiler::mFortProfs.end()) {  // error:  fname not found
    std::string estring("bl_proffortfuncstop error:  mFortProfs function not started:  ");
    estring += fName;
    if(Profiler::mFortProfsErrors.size() < mFortProfMaxErrors) {
      Profiler::mFortProfsErrors.push_back(estring);
    }
  } else {  // delete the pointer and remove fname from map
    delete it->second;
    Profiler::mFortProfs.erase(fName);
  }
}


#else

#include <BLFort.H>

BL_FORT_PROC_DECL(BL_PROFFORTFUNCSTART_CPP,bl_proffortfuncstart_cpp)
  (
   const int istr[], const int *NSTR
   )
{
}

BL_FORT_PROC_DECL(BL_PROFFORTFUNCSTOP_CPP,bl_proffortfuncstop_cpp)
  (
   const int istr[], const int *NSTR
   )
{
}


#endif





