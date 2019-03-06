#pragma once
#include <FEBioMech/FEMechModel.h>
#include "FECore/Timer.h"
#include "FECore/DataStore.h"
#include "FEBioPlot/PlotFile.h"
#include <FECore/FECoreKernel.h>
#include "febiolib_api.h"
#include <FECore/Logfile.h>

//-----------------------------------------------------------------------------
// Dump level determines the times the restart file is written
enum FE_Dump_Level {
	FE_DUMP_NEVER,			// never write a dump file
	FE_DUMP_MAJOR_ITRS,		// create a dump file at the end of each converged time step
	FE_DUMP_STEP			// create a dump file at the end of an analysis step
};

//-----------------------------------------------------------------------------
//! The FEBio model specializes the FEModel class to implement FEBio specific
//! functionality.
//!
//! In addition it adds support for all I/O capabilities. 
//!
class FEBIOLIB_API FEBioModel : public FEMechModel
{
public:
	//! constructor
	FEBioModel();

	//! destructor
	~FEBioModel();

	//! Initializes data structures
	bool Init() override;

	//! Resets data structures
	bool Reset() override;

	//! Solves the problem
	bool Solve() override;

public: // --- I/O functions ---

	//! input data from file
	bool Input(const char* szfile);

	//! write to plot file
	void Write(unsigned int nwhen);

	//! Write log data
	void WriteLog(unsigned int nwhen);

	//! write data to log file
	void WriteData();

	//! dump data to archive for restart
	void DumpData();

	//! add to log 
	void Log(int ntag, const char* szmsg) override;

	// get the log file
	Logfile& GetLogFile() { return m_log; }

public:
	//! set the problem title
	void SetTitle(const char* sz);

	//! get the problem title
	const std::string& GetTitle() const;

public: //! --- serialization for restarts ---
	
	//! Write or read data from archive
	void Serialize(DumpStream& ar) override;

protected:
	// helper functions for serialization
	void SerializeIOData   (DumpStream& ar);
	void SerializeDataStore(DumpStream& ar);

	bool InitLogFile();
	bool InitPlotFile();

	//! get a list of domains that belong to a specific material
	void DomainListFromMaterial(vector<int>& lmat, vector<int>& ldom);

public: // --- I/O functions ---
	//! Add data record
	void AddDataRecord(DataRecord* pd);

	//! Get the plot file
	PlotFile* GetPlotFile();

	// set the i/o files
	void SetInputFilename(const char* szfile);
	void SetLogFilename  (const char* szfile);
	void SetPlotFilename (const char* szfile);
	void SetDumpFilename (const char* szfile);

	//! Get the I/O file names
	const char* GetInputFileName();
	const char* GetLogfileName  ();
	const char* GetPlotFileName ();
	const char* GetDumpFileName ();

	//! get the file title
	const char* GetFileTitle();

public:
	double GetEndTime() const;

public: // Timers

	//! Return the total timer
	Timer& GetSolveTimer();

	//! return number of seconds of time spent in linear solver
	int GetLinearSolverTime() const;

public:
	//! set the debug level
	void SetDebugFlag(bool b) { m_debug = b; }

	//! get the debug level
	bool GetDebugFlag() { return m_debug; }

	//! set the dump level (for cold restarts)
	void SetDumpLevel(int dumpLevel) { m_dumpLevel = dumpLevel; }

	//! get the dump level
	int GetDumpLevel() const { return m_dumpLevel; }

private:
	void print_parameter(FEParam& p, int level = 0);
	void print_parameter_list(FEParameterList& pl, int level = 0);
	void print_parameter_list(FECoreBase* pc, int level = 0);
	void echo_input();

private:
	Timer		m_SolveTime;	//!< timer to track total time to solve problem
	Timer		m_InputTime;	//!< timer to track time to read model
	Timer		m_InitTime;		//!< timer to track model initialization
	Timer		m_IOTimer;		//!< timer to track output (include plot, dump, and data)

	PlotFile*	m_plot;			//!< the plot file
	bool		m_becho;		//!< echo input to logfile
	bool		m_debug;		//!< debug flag

	int			m_logLevel;		//!< output level for log file

	int			m_dumpLevel;	//!< level or writing restart file

private:
	// accumulative statistics
	int		m_ntimeSteps;		//!< total nr of time steps
	int		m_ntotalIters;		//!< total nr of equilibrium iterations
	int		m_ntotalRHS;		//!< total nr of right hand side evaluations
	int		m_ntotalReforms;	//!< total nr of stiffness reformations

protected: // file names
	char*	m_szfile_title;			//!< master input file title 
	char	m_szfile[MAX_STRING];	//!< master input file name (= path + title)
	char	m_szplot[MAX_STRING];	//!< plot output file name
	char	m_szlog [MAX_STRING];	//!< log output file name
	char	m_szdump[MAX_STRING];	//!< dump file name

	std::string	m_title;	//!< model title

private:
	Logfile	m_log;

	DECLARE_FECORE_CLASS();
};
