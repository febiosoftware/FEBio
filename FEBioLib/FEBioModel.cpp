#include "stdafx.h"
#include "FEBioModel.h"
#include "FEBioPlot/FEBioPlotFile.h"
#include "FEBioXML/FEBioImport.h"
#include "FEBioXML/FERestartImport.h"
#include "FECore/NodeDataRecord.h"
#include "FECore/ElementDataRecord.h"
#include <FEBioMech/ObjectDataRecord.h>
#include "FECore/NLConstraintDataRecord.h"
#include "FECore/log.h"
#include "FECore/FECoreKernel.h"
#include "FECore/DumpFile.h"
#include "FECore/DOFS.h"
#include <FECore/FEAnalysis.h>
#include <NumCore/MatrixTools.h>
#include <FECore/LinearSolver.h>
#include <FECore/FEDomain.h>
#include <FECore/FEMaterial.h>
#include "febio.h"
#include "version.h"
#include <iostream>
#include <fstream>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEBioModel, FEMechModel)
	ADD_PARAMETER(m_title   , "title"    );
	ADD_PARAMETER(m_logLevel, "log_level");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// echo the input data to the log file
extern void echo_input(FEBioModel& fem);

//-----------------------------------------------------------------------------
// Callback that guides FEBio output
bool output_cb(FEModel* pfem, unsigned int nwhen, void* pd)
{
	FEBioModel* pfebio = (FEBioModel*) pd;

	// write output to screen
	pfebio->WriteLog(nwhen);

	// write plot file
	pfebio->Write(nwhen);

	return true;
}

//-----------------------------------------------------------------------------
// Constructor of FEBioModel class.
FEBioModel::FEBioModel()
{
	m_logLevel = 1;

	m_dumpLevel = FE_DUMP_NEVER;

	// --- I/O-Data ---
	m_szfile_title = 0;
	m_szfile[0] = 0;
	m_szplot[0] = 0;
	m_szlog[0] = 0;
	m_szdump[0] = 0;
	m_debug = false;
	m_becho = true;
	m_plot = nullptr;
	m_writeMesh = false;

	m_ntimeSteps = 0;
	m_ntotalIters = 0;
	m_ntotalRHS = 0;
	m_ntotalReforms = 0;

	// Add the output callback
	// We call this function always since we want to flush the logfile for each event.
	AddCallback(output_cb, CB_ALWAYS, this);
}

//-----------------------------------------------------------------------------
FEBioModel::~FEBioModel()
{
	// close the plot file
	if (m_plot) { delete m_plot; m_plot = 0; }
	m_log.close();
}

//-----------------------------------------------------------------------------
Timer& FEBioModel::GetSolveTimer()
{
	return m_SolveTime;
}

//-----------------------------------------------------------------------------
//! return number of seconds of time spent in linear solver
int FEBioModel::GetLinearSolverTime()
{
	Timer* t = GetTimer(TimerID::Timer_Solve);
	return (int)t->peek();
}

//-----------------------------------------------------------------------------
//! set the debug level
void FEBioModel::SetDebugFlag(bool b) { m_debug = b; }

//! get the debug level
bool FEBioModel::GetDebugFlag() { return m_debug; }

//! set the dump level (for cold restarts)
void FEBioModel::SetDumpLevel(int dumpLevel) { m_dumpLevel = dumpLevel; }

//! get the dump level
int FEBioModel::GetDumpLevel() const { return m_dumpLevel; }

//-----------------------------------------------------------------------------
//! Set the title of the model
void FEBioModel::SetTitle(const char* sz)
{
	m_title = sz;
}

//-----------------------------------------------------------------------------
//! Return the title of the model
const std::string& FEBioModel::GetTitle() const
{
	return m_title;
}

//-----------------------------------------------------------------------------
double FEBioModel::GetEndTime() const
{
	return GetCurrentStep()->m_tend;
}

//=============================================================================
//
//		FEBioModel: I-O Functions
//
//=============================================================================

//-----------------------------------------------------------------------------
//! Add a data record to the data store
void FEBioModel::AddDataRecord(DataRecord* pd)
{
	DataStore& dataStore = GetDataStore();
	dataStore.AddRecord(pd); 
}

//-----------------------------------------------------------------------------
//! Get the plot file
PlotFile* FEBioModel::GetPlotFile()
{
	return m_plot;
}

//-----------------------------------------------------------------------------
//! Sets the name of the FEBio input file
void FEBioModel::SetInputFilename(const char* szfile)
{ 
	strcpy(m_szfile, szfile); 
	m_szfile_title = strrchr(m_szfile, '/');
	if (m_szfile_title == 0) 
	{
		m_szfile_title = strrchr(m_szfile, '\\'); 
		if (m_szfile_title == 0) m_szfile_title = m_szfile; else ++m_szfile_title;
	}
	else ++m_szfile_title;
}

//-----------------------------------------------------------------------------
//! Set the name of the log file
void FEBioModel::SetLogFilename(const char* szfile) 
{ 
	strcpy(m_szlog , szfile); 
}

//-----------------------------------------------------------------------------
//! Set the name of the plot file
void FEBioModel::SetPlotFilename(const char* szfile) 
{ 
	strcpy(m_szplot, szfile); 
}

//-----------------------------------------------------------------------------
//! Set the name of the restart archive (i.e. the dump file)
void FEBioModel::SetDumpFilename (const char* szfile) 
{ 
	strcpy(m_szdump, szfile); 
}

//-----------------------------------------------------------------------------
//! Return the name of the input file
const char* FEBioModel::GetInputFileName()
{ 
	return m_szfile; 
}

//-----------------------------------------------------------------------------
//! Return the name of the log file
const char* FEBioModel::GetLogfileName()
{ 
	return m_szlog;  
}

//-----------------------------------------------------------------------------
//! Return the name of the plot file
const char* FEBioModel::GetPlotFileName()
{ 
	return m_szplot; 
}

//-----------------------------------------------------------------------------
//! Return the dump file name.
const char* FEBioModel::GetDumpFileName()
{
	return	m_szdump;
}

//-----------------------------------------------------------------------------
//! get the file title (i.e. name of input file without the path)
const char* FEBioModel::GetFileTitle()
{ 
	return m_szfile_title; 
}

//=============================================================================
//    I N P U T
//=============================================================================

//-----------------------------------------------------------------------------
//! This routine reads in an input file and performs some initialization stuff.
//! The rest of the initialization is done in Init

bool FEBioModel::Input(const char* szfile)
{
	// start the timer
	TimerTracker t(&m_InputTime);

	// create file reader
	FEBioImport fim;

	feLog("Reading file %s ...", szfile);

	// Load the file
	if (fim.Load(*this, szfile) == false)
	{
		feLog("FAILED!\n");
		char szerr[256];
		fim.GetErrorMessage(szerr);
		feLog(szerr);

		return false;
	}
	else feLog("SUCCESS!\n");

	// set the input file name
	SetInputFilename(szfile);

	// see if user redefined output filenames
	if (fim.m_szdmp[0]) SetDumpFilename(fim.m_szdmp);
	if (fim.m_szlog[0]) SetLogFilename (fim.m_szlog);
	if (fim.m_szplt[0]) SetPlotFilename(fim.m_szplt);

	// set the plot file
	if (strcmp(fim.m_szplot_type, "febio") == 0)
	{
		FEBioPlotFile* pplt = new FEBioPlotFile(*this);
		m_plot = pplt;

		// set compression
		pplt->SetCompression(fim.m_nplot_compression);

		// define the plot file variables
		FEModel& fem = *GetFEModel();
		int NP = (int) fim.m_plot.size();
		for (int i=0; i<NP; ++i)
		{
			FEBioImport::PlotVariable& var = fim.m_plot[i];

			vector<int> item = var.m_item;
			if (item.empty() == false)
			{
				// TODO: currently, this is only supported for domain variables, where
				//       the list is a list of materials
				vector<int> lmat = var.m_item;

				// convert the material list to a domain list
				DomainListFromMaterial(lmat, item);
			}

			// add the plot output variable
			if (pplt->AddVariable(var.m_szvar, item, var.m_szdom) == false)
			{
				feLog("FATAL ERROR: Output variable \"%s\" is not defined\n", var.m_szvar);
				return false;
			}
		}
	}

	// add the data records
	int ND = (int)fim.m_data.size();
	for (int i=0; i<ND; ++i) AddDataRecord(fim.m_data[i]);

	// we're done reading
	return true;
}

//-----------------------------------------------------------------------------
//! This function finds all the domains that have a certain material
void FEBioModel::DomainListFromMaterial(vector<int>& lmat, vector<int>& ldom)
{
	FEMesh& mesh = GetMesh();

	// make sure the list is empty
	if (ldom.empty() == false) ldom.clear();

	// loop over all domains
	int ND = mesh.Domains();
	int NM = (int)lmat.size();
	for (int i = 0; i<ND; ++i)
	{
		FEDomain& di = mesh.Domain(i);
		int dmat = di.GetMaterial()->GetID();
		for (int j = 0; j<NM; ++j)
		{
			if (dmat == lmat[j])
			{
				ldom.push_back(i);
				break;
			}
		}
	}
}

//=============================================================================
//    O U T P U T
//=============================================================================

//! Write log data
void FEBioModel::WriteLog(unsigned int nwhen)
{
	FEAnalysis* step = GetCurrentStep();
	if (step == nullptr) return;

	if (nwhen == CB_STEP_SOLVED)
	{
		// output report
		feLog("\n\n N O N L I N E A R   I T E R A T I O N   I N F O R M A T I O N\n\n");
		feLog("\tNumber of time steps completed .................... : %d\n\n", step->m_ntimesteps);
		feLog("\tTotal number of equilibrium iterations ............ : %d\n\n", step->m_ntotiter);
		feLog("\tAverage number of equilibrium iterations .......... : %lg\n\n", (double)step->m_ntotiter / (double)step->m_ntimesteps);
		feLog("\tTotal number of right hand evaluations ............ : %d\n\n", step->m_ntotrhs);
		feLog("\tTotal number of stiffness reformations ............ : %d\n\n", step->m_ntotref);

		// print linear solver stats
		LinearSolver* ls = step->GetFESolver()->GetLinearSolver();
		if (ls)
		{
			LinearSolverStats stats = ls->GetStats();
			int nsolves = stats.backsolves;
			int niters = stats.iterations;
			double avgiters = (nsolves != 0 ? (double)niters / (double)nsolves : (double)niters);
			feLog("\n\n L I N E A R   S O L V E R   S T A T S\n\n");
			feLog("\tTotal calls to linear solver ........ : %d\n\n", nsolves);
			feLog("\tAvg iterations per solve ............ : %lg\n\n", avgiters);
		}

		// add to stats
		m_ntimeSteps += step->m_ntimesteps;
		m_ntotalIters += step->m_ntotiter;
		m_ntotalRHS += step->m_ntotrhs;
		m_ntotalReforms += step->m_ntotref;
	}

	if (nwhen == CB_SOLVED)
	{
		// for multistep analysis we'll print a grand total
		if (Steps() > 1)
		{
			feLog("\n\n N O N L I N E A R   I T E R A T I O N   S U M M A R Y\n\n");
			feLog("\tNumber of time steps completed .................... : %d\n\n", m_ntimeSteps);
			feLog("\tTotal number of equilibrium iterations ............ : %d\n\n", m_ntotalIters);
			feLog("\tTotal number of right hand evaluations ............ : %d\n\n", m_ntotalRHS);
			feLog("\tTotal number of stiffness reformations ............ : %d\n\n", m_ntotalReforms);
		}

		// get and print elapsed time
		char sztime[64];

		Timer* solveTimer = GetTimer(TimerID::Timer_Solve);
		solveTimer->time_str(sztime);
		feLog("\tTime in linear solver: %s\n\n", sztime);
	}

	// always flush the log
	m_log.flush();
}

//-----------------------------------------------------------------------------
//! Export state to plot file.
void FEBioModel::Write(unsigned int nwhen)
{
	TimerTracker t(&m_IOTimer);

	// get the current step
	FEAnalysis* pstep = GetCurrentStep();

	if (m_plot)
	{
		// get the plot level
		int nplt = pstep->GetPlotLevel();

		// if we don't want to plot anything we return
		if (nplt != FE_PLOT_NEVER)
		{
			// try to open the plot file
			if (nwhen == CB_STEP_ACTIVE)
			{
				if (m_plot->IsValid() == false)
				{
					if (m_plot->Open(*this, m_szplot) == false)
					{
						feLog("ERROR : Failed creating PLOT database\n");
						delete m_plot;
						m_plot = 0;
					}

					// Since it is assumed that for the first timestep
					// there are no loads or initial displacements, the case n=0 is skipped.
					// Therefor we can output those results here.
					// TODO: Offcourse we should actually check if this is indeed
					//       the case, otherwise we should also solve for t=0
					// Only output the initial state if requested
					if (m_plot)
					{
						bool bout = true;

						// if we're using the fixed time stepper, we check the plot range and zero state flag
						if (pstep->m_bautostep == false) bout = (pstep->m_nplotRange[0] == 0) || (pstep->m_bplotZero);

						// store initial time step (i.e. time step zero)
						double time = GetTime().currentTime;
						if (bout) m_plot->Write(*this, (float) time);
					}
				}
			}
			else
			{
				// assume we won't be writing anything
				bool bout = false;

				// see if we need to output something
				bool bdebug = GetDebugFlag();

				// write a new mesh section if needed
				if (nwhen == CB_REMESH)
				{
					m_writeMesh = true;
				}

				// when debugging we always output
				// (this could mean we may end up writing the same state multiple times)
				if (bdebug) bout = true;
				else
				{
					int currentStep = pstep->m_ntimesteps;
					int lastStep = pstep->m_ntime;
					int nmin = pstep->m_nplotRange[0]; if (nmin < 0) nmin = lastStep + nmin + 1;
					int nmax = pstep->m_nplotRange[1]; if (nmax < 0) nmax = lastStep + nmax + 1;

					bool inRange = true;
					bool isStride = true;
					if (pstep->m_bautostep == false)
					{
						inRange = false;
						if ((currentStep >= nmin) && (currentStep <= nmax)) inRange = true;

					}
                    isStride = ((pstep->m_ntimesteps - nmin) % pstep->m_nplot_stride) == 0;

					switch (nwhen)
					{
					case CB_MINOR_ITERS: if (nplt == FE_PLOT_MINOR_ITRS   ) bout = true; break;
					case CB_MAJOR_ITERS  : 
						if ((nplt == FE_PLOT_MAJOR_ITRS ) && inRange && isStride) bout = true; 
						if ((nplt == FE_PLOT_MUST_POINTS) && (pstep->m_timeController.m_nmust >= 0)) bout = true;
						if (nplt == FE_PLOT_AUGMENTATIONS) bout = true;
						break;
					case CB_AUGMENT: 
						if (nplt == FE_PLOT_AUGMENTATIONS) 
						{
							// Note that this is called before the augmentations.
							// The reason we store the state prior to the augmentations
							// is because the augmentations are going to change things such that
							// the system no longer in equilibrium. Since the model has to be converged
							// before we do augmentations, storing the model now will store an actual converged state.
							bout = true;
						}
						break;
					case CB_SOLVED : 
						if (nplt == FE_PLOT_FINAL) bout = true; 
						if (nplt == FE_PLOT_MAJOR_ITRS)
						{
							// we want to force storing the final time step, but we have to make sure 
							// it hasn't been stored already during the CB_MAJOR_ITERS callback.
							if ((inRange == false) || (isStride == false)) bout = true;
						}
						break;
					case CB_STEP_SOLVED: if (nplt == FE_PLOT_STEP_FINAL) bout = true;  break;
					}
				}

				// output the state if requested
				if (bout) 
				{
					// see if we need to write a new mesh section
					if (m_writeMesh) {
						FEBioPlotFile* plt = dynamic_cast<FEBioPlotFile*>(m_plot);
						plt->WriteMeshSection(*this);
						m_writeMesh = false;
					}

					double time = GetTime().currentTime;
					if (m_plot) m_plot->Write(*this, (float)time);
				}
			}
		}
	}

	// Dump converged state to the archive
	int ndump = GetDumpLevel();
	if (ndump != FE_DUMP_NEVER)
	{
		bool bdump = false;
		if ((nwhen == CB_STEP_SOLVED) && (ndump == FE_DUMP_STEP      )) bdump = true;
		if ((nwhen == CB_MAJOR_ITERS) && (ndump == FE_DUMP_MAJOR_ITRS)) bdump = true;
		if (bdump) DumpData();
	}

	// write the output data
	int nout = pstep->GetOutputLevel();
	if (nout != FE_OUTPUT_NEVER)
	{
		bool bout = false;
		switch (nwhen)
		{
		case CB_MINOR_ITERS: if (nout == FE_OUTPUT_MINOR_ITRS) bout = true; break;
		case CB_MAJOR_ITERS:
			if (nout == FE_OUTPUT_MAJOR_ITRS) bout = true;
			if ((nout == FE_OUTPUT_MUST_POINTS) && (pstep->m_timeController.m_nmust >= 0)) bout = true;
			break;
		case CB_SOLVED:
			if (nout == FE_OUTPUT_FINAL) bout = true;
			break;
		}

		if (bout) WriteData();
	}
}

//-----------------------------------------------------------------------------
//! Write user data to the logfile
void FEBioModel::WriteData()
{
	DataStore& dataStore = GetDataStore();
	dataStore.Write();
}

//-----------------------------------------------------------------------------
//! Dump state to archive for restarts
void FEBioModel::DumpData()
{
	DumpFile ar(*this);
	if (ar.Create(m_szdump) == false)
	{
		feLog("WARNING: Failed creating restart file (%s).\n", m_szdump);
	}
	else 
	{
		Serialize(ar);
		feLog("\nRestart point created. Archive name is %s\n", m_szdump);
	}
}

//-----------------------------------------------------------------------------
void FEBioModel::Log(int ntag, const char* szmsg)
{
	if      (ntag == 0) m_log.printf(szmsg);
	else if (ntag == 1) m_log.printbox("WARNING", szmsg);
	else if (ntag == 2) m_log.printbox("ERROR", szmsg);

	// Flushing the logfile each time we get here might be a bit overkill.
	// For now, I'm flushing the log file in the output_cb method.
//	m_log.flush();
}

//=============================================================================
//    R E S T A R T
//=============================================================================

class restart_exception : public std::runtime_error
{
public:
	restart_exception() : std::runtime_error("restart error") {}
	restart_exception(const char* msg) : std::runtime_error(msg) {}
};

//-----------------------------------------------------------------------------
//!  Reads or writes the current state to/from a binary file
//!  This is used to restart the solution from a saved position
//!  or to create a restart point.
//!  A version number is written to file to make sure the same
//!  format is used for reading and writing.
//! \param[in] ar the archive to which the data is serialized
//! \sa DumpFile
void FEBioModel::Serialize(DumpStream& ar)
{
	// don't need to do anything for running restarts
	if (ar.IsShallow())
	{
		// serialize model data
		FEMechModel::Serialize(ar);
	}
	else
	{
		if (ar.IsSaving())
		{
			// --- version number ---
			ar << (int) RSTRTVERSION;
		}
		else
		{
			// --- version ---
			int nversion;
			ar >> nversion;

			// make sure it is the right version
			if (nversion != RSTRTVERSION) throw restart_exception("incorrect version number");
		}

		// serialize model data
		FEMechModel::Serialize(ar);

		// serialize data store
		SerializeDataStore(ar);

		// --- Save IO Data
		SerializeIOData(ar);
	}
}

//-----------------------------------------------------------------------------
//! Serialization of FEBioModel data
void FEBioModel::SerializeIOData(DumpStream &ar)
{
	if (ar.IsSaving())
	{
		// file names
		ar << m_szfile << m_szplot << m_szlog << m_szdump;

		// plot file
		int npltfmt = 2;
		ar << npltfmt;

		// data records
		SerializeDataStore(ar);
	}
	else
	{
		// file names
		ar >> m_szfile >> m_szplot >> m_szlog >> m_szdump;

		// don't forget to call store the input file name so
		// that m_szfile_title gets initialized
		SetInputFilename(m_szfile);

		// get the plot file format (should be 2)
		int npltfmt = 0;
		ar >> npltfmt;
		assert(npltfmt == 2);

		// remove the plot file (if any)
		if (m_plot) { delete m_plot; m_plot = 0; }

		// create the plot file and open it for appending
		m_plot = new FEBioPlotFile(*this);
		if (m_plot->Append(*this, m_szplot) == false)
		{
			printf("FATAL ERROR: Failed reopening plot database %s\n", m_szplot);
			throw "FATAL ERROR";
		}

		// data records
		SerializeDataStore(ar);
	}
}

//-----------------------------------------------------------------------------
void FEBioModel::SerializeDataStore(DumpStream& ar)
{
	DataStore& dataStore = GetDataStore();
	if (ar.IsSaving())
	{
		int N = dataStore.Size();
		ar << N;
		for (int i=0; i<N; ++i)
		{
			DataRecord* pd = dataStore.GetDataRecord(i);

			int ntype = pd->m_type;
			ar << ntype;
			pd->Serialize(ar);
		}
	}
	else
	{
		int N;
		dataStore.Clear();
		ar >> N;
		for (int i=0; i<N; ++i)
		{
			int ntype;
			ar >> ntype;

			DataRecord* pd = 0;
			switch(ntype)
			{
			case FE_DATA_NODE: pd = new NodeDataRecord        (this, 0); break;
			case FE_DATA_ELEM: pd = new ElementDataRecord     (this, 0); break;
			case FE_DATA_RB  : pd = new ObjectDataRecord      (this, 0); break;
			case FE_DATA_NLC : pd = new NLConstraintDataRecord(this, 0); break;
			}
			assert(pd);
			pd->Serialize(ar);
			dataStore.AddRecord(pd);
		}
	}
}

//=============================================================================
//    I N I T I A L I Z A T I O N
//=============================================================================

//-----------------------------------------------------------------------------
//! This function performs one-time-initialization stuff. All the different 
//! modules are initialized here as well. This routine also performs some
//! data checks

bool FEBioModel::Init()
{
	TimerTracker t(&m_InitTime);

	// Open the logfile
	if (m_logLevel != 0)
	{
		if (InitLogFile() == false) return false;
	}

	// open plot database file
	FEAnalysis* step = GetCurrentStep();
	if (step->GetPlotLevel() != FE_PLOT_NEVER)
	{
		if (m_plot == 0) 
		{
			m_plot = new FEBioPlotFile(*this);
		}

		// see if a valid plot file name is defined.
		const char* szplt = GetPlotFileName();
		if (szplt[0] == 0)
		{
			// if not, we take the input file name and set the extension to .xplt
			char sz[1024] = {0};
			strcpy(sz, GetInputFileName());
			char *ch = strrchr(sz, '.');
			if (ch) *ch = 0;
			strcat(sz, ".xplt");
			SetPlotFilename(sz);
		}
	}

	// initialize model data
	if (FEMechModel::Init() == false) 
	{
		feLogError("Model initialization failed");
		return false;
	}

	// see if a valid dump file name is defined.
	const char* szdmp = this->GetDumpFileName();
	if (szdmp[0] == 0)
	{
		// if not, we take the input file name and set the extension to .dmp
		char sz[1024] = {0};
		strcpy(sz, GetInputFileName());
		char *ch = strrchr(sz, '.');
		if (ch) *ch = 0;
		strcat(sz, ".dmp");
		SetDumpFilename(sz);
	}

	// initialize data records
	DataStore& dataStore = GetDataStore();
	for (int i=0; i<dataStore.Size(); ++i)
	{
		if (dataStore.GetDataRecord(i)->Initialize() == false) return false;
	}

	// echo fem data to the logfile
	// we do this here (and not e.g. directly after input)
	// since the data can be changed after input, which is the case,
	// for instance, in the parameter optimization module
	if (m_becho) 
	{
		Logfile::MODE old_mode = m_log.GetMode();

		// don't output when no output is requested
		if (old_mode != Logfile::LOG_NEVER)
		{
			// we only output this data to the log file and not the screen
			m_log.SetMode(Logfile::LOG_FILE);

			// write output
			echo_input();

			// reset log mode
			m_log.SetMode(old_mode);
		}
	}

	// Alright, all initialization is done, so let's get busy !
	return true;
}

//-----------------------------------------------------------------------------
// Opens the log file
bool FEBioModel::InitLogFile()
{
	// Only do this if the log file is not valid
	if (!m_log.is_valid())
	{
		// see if a valid log file name is defined.
		const char* szlog = GetLogfileName();
		if (szlog[0] == 0)
		{
			// if not, we take the input file name and set the extension to .log
			char sz[1024] = {0};
			strcpy(sz, GetInputFileName());
			char *ch = strrchr(sz, '.');
			if (ch) *ch = 0;
			strcat(sz, ".log");
			SetLogFilename(sz);
		}

		// create a log stream
		LogFileStream* fp = new LogFileStream;
		m_log.SetFileStream(fp);
		if (fp->open(m_szlog) == false)
		{
			feLogError("Failed creating log file");
			return false;
		}

		// make sure we have a step
		FEAnalysis* step = GetCurrentStep();
		if (step == 0)
		{
			feLogError("No step defined.");
			return false;
		}

		// print welcome message to file
		Logfile::MODE m = m_log.SetMode(Logfile::LOG_FILE);
		febio::Hello(*fp);
		m_log.SetMode(m);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function resets the FEM data so that a new run can be done.
//! This routine is called from the optimization routine.

bool FEBioModel::Reset()
{
	// Reset model data
	FEMechModel::Reset();

	// re-initialize the log file
	if (m_logLevel != 0)
	{
		if (InitLogFile() == false) return false;
	}

	// open plot database file
	FEAnalysis* step =  GetCurrentStep();
	if (step->GetPlotLevel() != FE_PLOT_NEVER)
	{
		if (m_plot == 0) 
		{
			m_plot = new FEBioPlotFile(*this);
		}

		if (m_plot->Open(*this, m_szplot) == false)
		{
			feLogError("Failed creating PLOT database.");
			return false;
		}
	}

	m_ntimeSteps = 0;
	m_ntotalIters = 0;
	m_ntotalRHS = 0;
	m_ntotalReforms = 0;

	// do the callback
	DoCallback(CB_INIT);

	// All data is reset successfully
	return true;
}

//=============================================================================
//                               S O L V E
//=============================================================================

//-----------------------------------------------------------------------------
//! This is the main solve method. This function loops over all analysis steps
//! and solves each one in turn. 
//! \sa FEAnalysis

bool FEBioModel::Solve()
{
	// start the total time tracker
	m_SolveTime.start();

	// solve the FE model
	bool bconv = FEMechModel::Solve();

	// stop total time tracker
	m_SolveTime.stop();

	// print the elapsed time
	char sztime[64];
	m_SolveTime.time_str(sztime);
	feLog("\n Elapsed time : %s\n\n", sztime);

	// print additional stats to the log file only
	if (m_log.GetMode() & Logfile::LOG_FILE)
	{
		// print more detailed timing info to the log file
		Logfile::MODE old_mode = m_log.SetMode(Logfile::LOG_FILE);

		// sum up all the times spend in the linear solvers
		double total_time = 0.0;
		double input_time   = m_InputTime.GetTime(); total_time += input_time;
		double init_time    = m_InitTime.GetTime (); total_time += init_time;
		double solve_time   = m_SolveTime.GetTime(); total_time += solve_time;
		double io_time      = m_IOTimer.GetTime  ();
		double total_linsol = 0.0;
		double total_reform = 0.0;
		double total_stiff  = 0.0;
		double total_rhs    = 0.0;
		double total_update = 0.0;
		double total_qn     = 0.0;
		int NS = Steps();
		for (int i = 0; i<NS; ++i)
		{
			FEAnalysis* pstep = GetStep(i);
			FESolver* psolve = pstep->GetFESolver();
			if (psolve) 
			{
				total_linsol += GetTimer(TimerID::Timer_Solve    )->GetTime();
				total_reform += GetTimer(TimerID::Timer_Reform   )->GetTime();
				total_stiff  += GetTimer(TimerID::Timer_Stiffness)->GetTime();
				total_rhs    += GetTimer(TimerID::Timer_Residual )->GetTime();
				total_update += GetTimer(TimerID::Timer_Update   )->GetTime();
				total_qn     += GetTimer(TimerID::Timer_QNUpdate )->GetTime();
			}
		}

		feLog(" T I M I N G   I N F O R M A T I O N\n\n");
		Timer::time_str(input_time  , sztime); feLog("\tInput time ...................... : %s (%lg sec)\n\n", sztime, input_time  );
		Timer::time_str(init_time   , sztime); feLog("\tInitialization time ............. : %s (%lg sec)\n\n", sztime, init_time   );
		Timer::time_str(solve_time  , sztime); feLog("\tSolve time ...................... : %s (%lg sec)\n\n", sztime, solve_time  );
		Timer::time_str(io_time     , sztime); feLog("\t   IO-time (plot, dmp, data) .... : %s (%lg sec)\n\n", sztime, io_time     );
		Timer::time_str(total_reform, sztime); feLog("\t   reforming stiffness .......... : %s (%lg sec)\n\n", sztime, total_reform);
		Timer::time_str(total_stiff , sztime); feLog("\t   evaluating stiffness ......... : %s (%lg sec)\n\n", sztime, total_stiff );
		Timer::time_str(total_rhs   , sztime); feLog("\t   evaluating residual .......... : %s (%lg sec)\n\n", sztime, total_rhs   );
		Timer::time_str(total_update, sztime); feLog("\t   model update ................. : %s (%lg sec)\n\n", sztime, total_update);
		Timer::time_str(total_qn    , sztime); feLog("\t   QN updates ................... : %s (%lg sec)\n\n", sztime, total_qn);
		Timer::time_str(total_linsol, sztime); feLog("\t   time in linear solver ........ : %s (%lg sec)\n\n", sztime, total_linsol);
		Timer::time_str(total_time  , sztime); feLog("\tTotal elapsed time .............. : %s (%lg sec)\n\n", sztime, total_time  );


		m_log.SetMode(old_mode);

		if (bconv)
		{
			feLog("\n N O R M A L   T E R M I N A T I O N\n\n");
		}
		else
		{
			feLog("\n E R R O R   T E R M I N A T I O N\n\n");
		}

		// flush the log file
		m_log.flush();
	}

	// close the plot file
	if (m_plot) m_plot->Close();

	// We're done !
	return bconv;
}
