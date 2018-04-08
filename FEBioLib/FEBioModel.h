#pragma once
#include "FECore/FEModel.h"
#include "FECore/Timer.h"
#include "FECore/DataStore.h"
#include "FEBioPlot/PlotFile.h"
#include <FECore/FECoreKernel.h>

//-----------------------------------------------------------------------------
//! The FEBio model specializes the FEModel class to implement FEBio specific
//! functionality.
//!
//! In addition it adds support for all I/O capabilities. 
//!
class FECOREDLL_EXPORT FEBioModel : public FEModel
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

public:
	//! set the problem title
	void SetTitle(const char* sz);

	//! get the problem title
	const char* GetTitle();

public: //! --- serialization for restarts ---
	
	//! Write or read data from archive
	void Serialize(DumpStream& ar) override;

protected:
	// helper functions for serialization
	void SerializeIOData   (DumpStream& ar);
	void SerializeDataStore(DumpStream& ar);

	bool InitLogFile();
	bool InitPlotFile();

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

	//! return the data store
	DataStore& GetDataStore();

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

private:
	Timer		m_SolveTime;	//!< timer to track total time to solve problem
	Timer		m_InputTime;	//!< timer to track time to read model
	Timer		m_InitTime;		//!< timer to track model initialization
	Timer		m_IOTimer;		//!< timer to track output (include plot, dump, and data)

	DataStore	m_Data;			//!< the data store used for data logging
	PlotFile*	m_plot;			//!< the plot file
	bool		m_becho;		//!< echo input to logfile \todo Make this a command line option
	bool		m_debug;		//!< debug flag

protected: // file names
	char*	m_szfile_title;			//!< master input file title 
	char	m_szfile[MAX_STRING];	//!< master input file name (= path + title)
	char	m_szplot[MAX_STRING];	//!< plot output file name
	char	m_szlog [MAX_STRING];	//!< log output file name
	char	m_szdump[MAX_STRING];	//!< dump file name

	char	m_sztitle[MAX_STRING];	//!< model title

	DECLARE_PARAMETER_LIST();
};
