#pragma once
#include "FECore/FEModel.h"
#include "DataStore.h"
#include "Timer.h"
#include "FEElasticMixture.h"
#include "FEUncoupledElasticMixture.h"

//-----------------------------------------------------------------------------
//!< The FEBio model specializes the FEModel class to implement FEBio specific
//!< functionality.
class FEBioModel : public FEModel
{
public:
	//!< constructor
	FEBioModel();

	//! input data from file
	bool Input(const char* szfile);

	//! Initializes data structures
	bool Init();

	//! Resets data structures
	bool Reset();

	//! Solves the problem
	bool Solve();

	//! find a boundary condition from the ID
	FEBoundaryCondition* FindBC(int nid);

public: // --- virtual I/O functions ---
	//! write to plot file
	void Write();

	//! write data to log file
	void WriteData();

	//! dump data to archive for restart
	void DumpData();

protected: //! --- serialization for restarts ---
	
	//! Write or read data from archive
	bool Serialize(DumpFile& ar);

	// helper functions for serialization
	void SerializeMaterials   (DumpFile& ar);
	void SerializeAnalysisData(DumpFile& ar);
	void SerializeGeometry    (DumpFile& ar);
	void SerializeMesh        (DumpFile& ar);
	void SerializeContactData (DumpFile& ar);
	void SerializeBoundaryData(DumpFile& ar);
	void SerializeIOData      (DumpFile& ar);
	void SerializeLoadData    (DumpFile& ar);
	void SerializeConstants   (DumpFile& ar);
	void SerializeDataStore   (DumpFile& ar);

public: //! --- parameter functions ---

	//! return a pointer to the named variable
	double* FindParameter(const char* szname);

	//! Evaluate parameter list
	void EvaluateParameterList(FEParameterList& pl);
	void EvaluateMaterialParameters(FEMaterial* pm);

protected: // --- initialization functions ---

	//! Initialize mesh data
	bool InitMesh();

	//! Initialize rigid bodies
	bool CreateRigidBodies();

	//! Initialize material data
	bool InitMaterials();

	//! Initializes contact data
	bool InitContact();

	//! Initialize poroelastic/biphasic and solute data
	bool InitPoroSolute();

public: // --- I/O functions ---
	//! Add data record
	void AddDataRecord(DataRecord* pd) { m_Data.AddRecord(pd); }

	//! Set the plot file
	void SetPlotFile(PlotFile* pplt) { m_plot = pplt; }

	// set the i/o files
	void SetInputFilename(const char* szfile);
	void SetLogFilename  (const char* szfile) { strcpy(m_szlog , szfile); }
	void SetPlotFilename (const char* szfile) { strcpy(m_szplot, szfile); }
	void SetDumpFilename (const char* szfile) { strcpy(m_szdump, szfile); }

	void SetPlotFileNameExtension(const char* szext);

	//! Get the I/O file names
	const char* GetInputFileName() { return m_szfile; }
	const char* GetLogfileName  () { return m_szlog;  }
	const char* GetPlotFileName () { return m_szplot; }

	//! get the file title
	const char* GetFileTitle() { return m_szfile_title; }

public:
	Timer		m_TotalTime;	//!< Create timer to track total running time
	PlotFile*	m_plot;			//!< the plot file
	DataStore	m_Data;			//!< the data store used for data logging
	bool		m_becho;		//!< echo input to logfile (TODO: Make this a command line option)

protected:
	// file names
	char*	m_szfile_title;			//!< master input file title 
	char	m_szfile[MAX_STRING];	//!< master input file name (= path + title)
	char	m_szplot[MAX_STRING];	//!< plot output file name
	char	m_szlog [MAX_STRING];	//!< log output file name
	char	m_szdump[MAX_STRING];	//!< dump file name
};
