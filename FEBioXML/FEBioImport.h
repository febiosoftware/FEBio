#pragma once

#include "FileImport.h"
#include "XMLReader.h"
#include "FEBioLib/FEAnalysisStep.h"
#include "FEBioLib/FESolver.h"
#include <map>
#include <string>
using namespace std;

class FEFEBioImport;

//-----------------------------------------------------------------------------
// Base class for XML sections parsers
class FEBioFileSection
{
public:
	FEBioFileSection(FEFEBioImport* pim) { m_pim = pim; }

	virtual void Parse(XMLTag& tag) = 0;

	FEModel* GetFEModel();
	FEAnalysisStep* GetStep();

protected:
	FEFEBioImport*	m_pim;
};

//-----------------------------------------------------------------------------
// class that manages file section parsers
class FEBioFileSectionMap : public map<string, FEBioFileSection*>
{
public:
	~FEBioFileSectionMap();
};

//=============================================================================
//! Implements a class to import FEBio input files
//!
class FEFEBioImport : public FEFileImport
{
public:
	// Element types
	enum { ET_HEX8, ET_HEX20, ET_PENTA6, ET_TET4, ET_UT4, ET_TETG1, ET_TET10, ET_QUAD4, ET_TRI3, ET_TRUSS2 };

	// element classes
	enum { EC_STRUCT, EC_RIGID, EC_PORO, EC_HEAT };

public:
	class InvalidVersion{};
	class InvalidMaterial
	{ 
	public: 
		InvalidMaterial(int nel) : m_nel(nel){}
		int m_nel; 
	};
	class InvalidDomainType{};
	class FailedCreatingDomain{};
	class InvalidElementType{};
	class FailedLoadingPlugin
	{
	public:
		FailedLoadingPlugin(const char* sz) : m_szfile(sz) {}
		const char* FileName() { return m_szfile.c_str(); }
	public:
		string	m_szfile;
	};
	class DuplicateMaterialSection {};

public:
	//! Load the model data from file.
	bool Load(FEModel& fem, const char* szfile);

public:
	FEModel* GetFEModel() { return m_pfem; }
	FEMesh* GetFEMesh() { return m_pMesh; }
	FEAnalysis*	GetStep();

	int Version() { return m_nversion; }

	bool ReadParameter(XMLTag& tag, FEParameterList& pl);

	void ReadList(XMLTag& tag, vector<int>& l);

	FEAnalysis* CreateNewStep();

public:
	void SetDumpfileName(const char* sz) { sprintf(m_szdmp, sz); }
	void SetLogfileName (const char* sz) { sprintf(m_szlog, sz); }
	void SetPlotfileName(const char* sz) { sprintf(m_szplt, sz); }

protected:
	void ParseVersion			(XMLTag& tag);

public:
	FEModel*	m_pfem;		//!< pointer to the fem class
	FEAnalysis*	m_pStep;	//!< pointer to current analysis step
	FEMesh*		m_pMesh;	//!< pointer to the mesh class

public:
	char	m_szdmp[256];
	char	m_szlog[256];
	char	m_szplt[256];

public:
	int	m_ntet4;	// tetrahedral integration rule
	int	m_nut4;		// integration rule for stabilization of UT4
	int m_nsteps;	// nr of step sections read
	int	m_nstep_type;	// step type
	int	m_nmat;		// nr of materials

	bool	m_b3field;	// three-field element flag
	int		m_nhex8;	// hex integration rule

protected:
	int	m_nversion;	// version of file
};
