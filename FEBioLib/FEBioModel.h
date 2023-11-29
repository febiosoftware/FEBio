/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include <FEBioMech/FEMechModel.h>
#include <FECore/Timer.h>
#include <FECore/DataStore.h>
#include <FEBioPlot/PlotFile.h>
#include <FECore/FECoreKernel.h>
#include <FEBioLib/Logfile.h>
#include "febiolib_api.h"

//-----------------------------------------------------------------------------
// Dump level determines the times the restart file is written
enum FE_Dump_Level {
	FE_DUMP_NEVER,			// never write a dump file
	FE_DUMP_MAJOR_ITRS,		// create a dump file at the end of each converged time step
	FE_DUMP_STEP,			// create a dump file at the end of an analysis step
	FE_DUMP_MUST_POINTS     // create a dump file only on must-points
};

//-----------------------------------------------------------------------------
struct ModelStats {
	int		ntimeSteps;		//!< total nr of time steps
	int		ntotalIters;	//!< total nr of equilibrium iterations
	int		ntotalRHS;		//!< total nr of right hand side evaluations
	int		ntotalReforms;	//!< total nr of stiffness reformations
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

public: // --- I/O functions ---

	//! input data from file
	bool Input(const char* szfile);

	//! handle output
	void Write(unsigned int nwhen);

	// write to plot file
	void WritePlot(unsigned int nevent);

	//! write data to log file
	void WriteData(unsigned int nevent);

	//! dump data to archive for restart
	void DumpData(int nevent);

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

private:
	static bool handleCB(FEModel* fem, unsigned int nwhen, void* pd);
	bool processEvent(int nevent);

	void on_cb_solved();
	void on_cb_stepSolved();

protected:
	// helper functions for serialization
	void SerializeIOData   (DumpStream& ar);
	void SerializeDataStore(DumpStream& ar);
	void SerializePlotData (DumpStream& ar);

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
	void SetInputFilename(const std::string& sfile);
	void SetLogFilename  (const std::string& sfile);
	void SetPlotFilename (const std::string& sfile);
	void SetDumpFilename (const std::string& sfile);

	//! Get the I/O file names
	const std::string& GetInputFileName();
	const std::string& GetLogfileName  ();
	const std::string& GetPlotFileName ();
	const std::string& GetDumpFileName ();

	//! get the file title
	const std::string& GetFileTitle();

	// set append-on-restart flag
	void SetAppendOnRestart(bool b);
	bool AppendOnRestart() const;

public:
	double GetEndTime() const;

public: // Timers

	//! Return the total timer
	Timer& GetSolveTimer();

	//! return number of seconds of time spent in linear solver
	int GetLinearSolverTime();

public:
	//! set the debug level
	void SetDebugLevel(int debugLvl);

	//! get the debug level
	int GetDebugLevel();

	//! set the dump level (for cold restarts)
	void SetDumpLevel(int dumpLevel);

	//! get the dump level
	int GetDumpLevel() const;

	//! Set the dump stride
	void SetDumpStride(int n);

	//! get the dump stride
	int GetDumpStride() const;

	//! Set the log level
	void SetLogLevel(int logLevel);

	//! Get the stats 
	ModelStats GetModelStats() const;

	// flag to show warnings and errors
	void ShowWarningsAndErrors(bool b);
	bool ShowWarningsAndErrors() const;

private:
	void print_parameter(FEParam& p, int level = 0);
	void print_parameter_list(FEParameterList& pl, int level = 0);
	void print_parameter_list(FECoreBase* pc, int level = 0);
	void echo_input();

private:
	void UpdatePlotObjects();

private:
	Timer		m_InputTime;	//!< timer to track time to read model
	Timer		m_InitTime;		//!< timer to track model initialization
	Timer		m_IOTimer;		//!< timer to track output (include plot, dump, and data)

	PlotFile*	m_plot;			//!< the plot file
	bool		m_becho;		//!< echo input to logfile
	int			m_ndebug;		//!< debug level flag
	bool		m_writeMesh;	//!< write a new mesh section

	bool		m_bshowErrors;	//!< print warnings and errors

	int			m_logLevel;		//!< output level for log file

	int			m_dumpLevel;	//!< level or writing restart file
	int			m_dumpStride;	//!< write dump file every nth iterations

private:
	// accumulative statistics
	ModelStats	m_stats;

protected: // file names
	std::string		m_sfile_title;		//!< input file title 
	std::string		m_sfile;			//!< input file name (= path + title)
	std::string		m_splot;			//!< plot output file name
	std::string		m_slog ;			//!< log output file name
	std::string		m_sdump;			//!< dump file name

	std::string	m_title;	//!< model title

protected:
	bool					m_pltAppendOnRestart;
	int						m_lastUpdate;

private:
	Logfile	m_log;

	DECLARE_FECORE_CLASS();
};
