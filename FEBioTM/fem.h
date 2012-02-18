#pragma once
#include "FEBioLib/FEBioModel.h"
#include "FEBioPlot/PlotFile.h"
#include "FEBioLib/Timer.h"
#include "FEBioLib/DataStore.h"

//-----------------------------------------------------------------------------
// The FE model class
class FEM : public FEBioModel
{
public:
	// constructor
	FEM();

	// Initialization
	virtual bool Init() { return false; }

	//! Resets data structures
	virtual bool Reset() { return false; }

	// solve the model
	virtual bool Solve() { return false; }

	virtual void PushState() {}
	virtual void PopState () {}

	//! Evaluate parameter list
	virtual void EvaluateParameterList(FEParameterList& pl) {}
	virtual void EvaluateMaterialParameters(FEMaterial* pm) {}

	//! return a pointer to the named variable
	virtual double* FindParameter(const char* szname) { return 0; }

	//! find a boundary condition from the ID
	virtual FEBoundaryCondition* FindBC(int nid) { return 0; }

	//! check for user interruption
	virtual void CheckInterruption() {}

	// input from file
	virtual bool Input(const char* szfile);

	//! write to plot file
	virtual void Write() {}

	//! write data to log file
	virtual void WriteData() {}

	//! write data to dump file
	virtual void DumpData() {}

	//! serialize data
	virtual bool Serialize(DumpFile& ar) { return false; }

	//! Add data record
	virtual void AddDataRecord(DataRecord* pd) {}

	//! set plot file
	void SetPlotFile(PlotFile* pplt) {}

	virtual void SetPlotFileNameExtension(const char *szext) {}

	// set the i/o files
	void SetInputFilename(const char* szfile);
	void SetLogFilename  (const char* szfile) { strcpy(m_szlog , szfile); }
	void SetPlotFilename (const char* szfile) { strcpy(m_szplot, szfile); }
	void SetDumpFilename (const char* szfile) { strcpy(m_szdump, szfile); }

public:
	Timer	m_TotalTime;		//!< Create timer to track total running time

public:
	// --- I/O-Data --- 
	//{
		PlotFile*	m_plot;		//!< the plot file
		DataStore	m_Data;		//!< the data store used for data logging
		
		bool		m_becho;			//!< echo input to logfile (TODO: Make this a command line option)

protected:
		// file names
		char*	m_szfile_title;			//!< master input file title 
		char	m_szfile[MAX_STRING];	//!< master input file name (= path + title)
		char	m_szplot[MAX_STRING];	//!< plot output file name
		char	m_szlog [MAX_STRING];	//!< log output file name
		char	m_szdump[MAX_STRING];	//!< dump file name
	//}
};
