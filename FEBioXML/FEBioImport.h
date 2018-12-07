#pragma once
#include "FileImport.h"
#include "XMLReader.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FESolver.h"
#include "FECore/DataStore.h"
#include <FECore/FEMesh.h>
#include <FECore/FESurfaceMap.h>
#include <FECore/tens3d.h>
#include <string>
using namespace std;

class FENodeSet;

class FEBioImport;

//-----------------------------------------------------------------------------
//! Base class for FEBio (feb) file sections.
class FEBioFileSection : public FEFileSection
{
public:
	FEBioFileSection(FEBioImport* feb);

	FEBioImport* GetFEBioImport();
};

//=============================================================================
//! Implements a class to import FEBio input files
//!
class FEBioImport : public FEFileImport
{
public:
	// invalid version
	class InvalidVersion : public FEFileException
	{
	public: InvalidVersion();
	};

	// invalid material defintion
	class InvalidMaterial : public FEFileException
	{ 
	public: InvalidMaterial(int nel);
	};

	// invalid domain type
	class InvalidDomainType : public FEFileException
	{
	public: InvalidDomainType();
	};

	// cannot create domain
	class FailedCreatingDomain : public FEFileException
	{
	public: FailedCreatingDomain();
	};

	// invalid element type
	class InvalidElementType : public FEFileException
	{
	public: InvalidElementType();
	};

	// failed loading plugin
	class FailedLoadingPlugin : public FEFileException
	{
	public: FailedLoadingPlugin(const char* sz);
	};

	// duplicate material section
	class DuplicateMaterialSection : public FEFileException
	{
	public: DuplicateMaterialSection();
	};

	// invalid domain material
	class InvalidDomainMaterial : public FEFileException
	{ 
	public: InvalidDomainMaterial();
	};

	// missing property
	class MissingProperty : public FEFileException
	{
	public: MissingProperty(const std::string& matName, const char* szprop);
	};

	//! Failed allocating solver
	class FailedAllocatingSolver : public FEFileException
	{
	public: FailedAllocatingSolver(const char* sztype);
	};

	//! error in data generation
	class DataGeneratorError : public FEFileException
	{
	public: 
		DataGeneratorError();
	};

	//! failed building a part
	class FailedBuildingPart : public FEFileException
	{
	public:
		FailedBuildingPart(const std::string& partName);
	};

	// Error while reading mesh data section
	class MeshDataError : public FEFileException
	{
	public:
		MeshDataError();
	};

public:
	//-------------------------------------------------------------------------
	class PlotVariable
	{
	public:
		PlotVariable(const PlotVariable& pv);
        PlotVariable(const std::string& var, vector<int>& item, const char* szdom = "");
        
	public:
		char		m_szvar[128];	//!< name of output variable
        char        m_szdom[128];    //!< (optional) name of domain
		vector<int>	m_item;			//!< (optional) list of items
	};

public:
	//! constructor
	FEBioImport();

	//! destructor
	~FEBioImport();

	//! open the file
	bool Load(FEModel& fem, const char* szfile);

	//! read the contents of a file
	bool ReadFile(const char* szfile, bool broot = true);

public:
	FEMesh* GetFEMesh() { return m_pMesh; }

public:
	void SetDumpfileName(const char* sz);
	void SetLogfileName (const char* sz);
	void SetPlotfileName(const char* sz);

    void AddPlotVariable(const char* szvar, vector<int>& item, const char* szdom = "");

	void SetPlotCompression(int n);
    
	void AddDataRecord(DataRecord* pd);

public:
	// Helper functions for reading node sets, surfaces, etc.
	FENodeSet* ParseNodeSet(XMLTag& tag, const char* szatt = "set");
	FESurface* ParseSurface(XMLTag& tag, const char* szatt = "surf");

protected:
	void ParseVersion(XMLTag& tag);

	void BuildFileSectionMap(int nversion);

public:
	FEMesh*			m_pMesh;	//!< pointer to the mesh class

public:
	char	m_szdmp[512];
	char	m_szlog[512];
	char	m_szplt[512];

public:
	char					m_szplot_type[256];
	vector<PlotVariable>	m_plot;
	int						m_nplot_compression;

	vector<DataRecord*>		m_data;
};
