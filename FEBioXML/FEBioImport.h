#pragma once
#include "FileImport.h"
#include "XMLReader.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FESolver.h"
#include "FECore/DataStore.h"
#include <map>
#include <string>
using namespace std;

class FEBioImport;

//-----------------------------------------------------------------------------
// Base class for XML sections parsers
class FEBioFileSection
{
public:
	FEBioFileSection(FEBioImport* pim) { m_pim = pim; }

	virtual void Parse(XMLTag& tag) = 0;

	FEModel* GetFEModel();
	FECore::FEAnalysis* GetStep();

protected:
	FEBioImport*	m_pim;
};

//-----------------------------------------------------------------------------
// class that manages file section parsers
class FEBioFileSectionMap : public map<string, FEBioFileSection*>
{
public:
	~FEBioFileSectionMap();
	void Clear();
};

//=============================================================================
//! Implements a class to import FEBio input files
//!
class FEBioImport : public FEFileImport
{
public:
	// Base class for FEBio import exceptions
	// Derived classes should set the error string in their constructor
	class Exception
	{
	public:
		enum {MAX_ERR_STRING = 1024};
	public:
		// retrieve the error string
		const char* GetErrorString() const { return m_szerr; }
	protected:
		// set the error string (used by derived classes)
		void SetErrorString(const char* sz, ...);
	protected:
		char	m_szerr[MAX_ERR_STRING];
	};

	// invalid version
	class InvalidVersion : public Exception
	{
	public: InvalidVersion();
	};

	// invalid material defintion
	class InvalidMaterial : public Exception
	{ 
	public: InvalidMaterial(int nel);
	};

	// invalid domain type
	class InvalidDomainType : public Exception
	{
	public: InvalidDomainType();
	};

	// cannot create domain
	class FailedCreatingDomain : public Exception
	{
	public: FailedCreatingDomain();
	};

	// invalid element type
	class InvalidElementType : public Exception
	{
	public: InvalidElementType();
	};

	// failed loading plugin
	class FailedLoadingPlugin : public Exception
	{
	public: FailedLoadingPlugin(const char* sz);
	};

	// duplicate material section
	class DuplicateMaterialSection : public Exception
	{
	public: DuplicateMaterialSection();
	};

	// invalid domain material
	class InvalidDomainMaterial : public Exception
	{ 
	public: InvalidDomainMaterial();
	};

	// missing material property
	class MissingMaterialProperty : public Exception
	{
	public: MissingMaterialProperty(const char* szmat, const char* szprop);
	};

	//! Failed allocating solver
	class FailedAllocatingSolver : public Exception
	{
	public: FailedAllocatingSolver(const char* sztype);
	};

	//! Invalid node ID
	class InvalidNodeID : public Exception
	{
	public: InvalidNodeID();
	};

	//! missing slave surface
	class MissingSlaveSurface : public Exception
	{
	public:
		MissingSlaveSurface();
	};

	//! missing master surface
	class MissingMasterSurface : public Exception
	{
	public:
		MissingMasterSurface();
	};

public:
	//-------------------------------------------------------------------------
	class PlotVariable
	{
	public:
		PlotVariable(const PlotVariable& pv);
        PlotVariable(const char* szvar, vector<int>& item, const char* szdom = "");
        
	public:
		char		m_szvar[64];	//!< name of output variable
        char        m_szdom[64];    //!< (optional) name of domain
		vector<int>	m_item;			//!< (optional) list of items
	};

	//-------------------------------------------------------------------------
	class XMLParam
	{
		enum {MAX_TAG = 128};
	public:
		XMLParam() { m_szname[0] = 0; m_szval[0] = 0; }
	public:
		char	m_szname[MAX_TAG];	// parameter name
		char	m_szval[MAX_TAG];	// parameter value
	};


public:
	//! constructor
	FEBioImport();

	//! Load the model data from file.
	bool Load(FEModel& fem, const char* szfile);

	//! read the contents of a file
	bool ReadFile(const char* szfile, bool broot = true);

public:
	FEModel* GetFEModel() { return m_pfem; }
	FEMesh* GetFEMesh() { return m_pMesh; }
	FECore::FEAnalysis*	GetStep();

	int Version() { return m_nversion; }

	bool ReadParameter(XMLTag& tag, FEParameterList& pl, const char* szparam = 0);
	bool ReadParameter(XMLTag& tag, FECoreBase* pc, const char* szparam = 0);

	void ReadList(XMLTag& tag, vector<int>& l);

	FESolver* BuildSolver(const char* sztype, FEModel& fem);
	FECore::FEAnalysis* CreateNewStep();

public:
	void SetDumpfileName(const char* sz);
	void SetLogfileName (const char* sz);
	void SetPlotfileName(const char* sz);

    void AddPlotVariable(const char* szvar, vector<int>& item, const char* szdom = "");

	void SetPlotCompression(int n);
    
	void AddDataRecord(DataRecord* pd);

public:
	// Helper functions for reading node sets, surfaces, etc.
	FENodeSet* ParseNodeSet(XMLTag& tag);

public:
	void ClearParams();

	XMLParam* FindParameter(const char* sz);

	void AddParameter(const char* szname, const char* szval);

public:
	const char* get_value_string(XMLTag& tag);
	void value(XMLTag& tag, int&    n);
	void value(XMLTag& tag, double& g);
	void value(XMLTag& tag, bool&   b);
	void value(XMLTag& tag, vec3d&  v);
	void value(XMLTag& tag, mat3d&  m);
	void value(XMLTag& tag, mat3ds& m);
	void value(XMLTag& tag, char* szstr);
	int value(XMLTag& tag, int* pi, int n);
	int value(XMLTag& tag, double* pf, int n);

protected:
	void ParseVersion			(XMLTag& tag);

public:
	FEModel*			m_pfem;		//!< pointer to the fem class
	FECore::FEAnalysis*	m_pStep;	//!< pointer to current analysis step
	FEMesh*				m_pMesh;	//!< pointer to the mesh class

public:
	char	m_szpath[512];
	char	m_szdmp[512];
	char	m_szlog[512];
	char	m_szplt[512];

public:
	char					m_szplot_type[256];
	vector<PlotVariable>	m_plot;
	int						m_nplot_compression;

	vector<DataRecord*>		m_data;

	vector<XMLParam>	m_Param;	// parameter list

public:
	char	m_szmod[256];	//!< module type string
	int		m_nsteps;		//!< nr of step sections read
	int		m_maxid;		//!< max element ID

	bool	m_b3field_hex;	//!< three-field element flag for hex (and wedge elements)
	bool	m_b3field_tet;	//!< three-field element flag for quadratic tets
	bool	m_but4;		//!< use UT4 formulation flag
	FE_Element_Type		m_nhex8;	//!< hex integration rule
	FE_Element_Type		m_ntet4;	//!< tet4 integration rule
	FE_Element_Type		m_ntet10;	//!< tet10 integration rule
	FE_Element_Type		m_ntet15;	//!< tet15 integration rule
	FE_Element_Type		m_ntri6;	//!< tri6 integration rule
	FE_Element_Type		m_ntri7;	//!< tri7 integration rule
	FE_Element_Type		m_ntri3;	//!< tri3 integration rule
    FE_Element_Type		m_nquad4;	//!< quad4 integration rule
    FE_Element_Type		m_nquad8;	//!< quad8 integration rule
    FE_Element_Type		m_nquad9;	//!< quad9 integration rule
    FE_Element_Type		m_nfquad4;	//!< fquad4 integration rule

public:
	// Build the node ID table
	void BuildNodeList();

	// find a node index from its ID
	int FindNodeFromID(int nid);

	// convert an array of nodal ID to nodal indices
	void GlobalToLocalID(int* l, int n, vector<int>& m);

	//! read a nodal ID
	//! This assumes the node ID is defined via the "id" attribute
	int ReadNodeID(XMLTag& tag);

protected:
	int			m_node_off;		//!< node offset (i.e. lowest node ID)
	vector<int>	m_node_list;	//!< map node ID's to their nodes.

protected:
	FEBioFileSectionMap	m_map;

protected:
	int	m_nversion;	// version of file
};
