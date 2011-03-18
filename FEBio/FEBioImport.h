#pragma once

#include "FileImport.h"
#include "XMLReader.h"
#include "FETransverselyIsotropic.h"
#include "FERigid.h"
#include "FEElasticMixture.h"
#include "FEUncoupledElasticMixture.h"
#include <map>
#include <string>
using namespace std;

class FEFEBioImport;

//-----------------------------------------------------------------------------
// Base class for XML sections parsers
class FileSection
{
public:
	FileSection(FEFEBioImport* pim) { m_pim = pim; }

	virtual void Parse(XMLTag& tag) = 0;

	FEM* GetFEM();
	FEAnalysis* GetStep();

protected:
	FEFEBioImport*	m_pim;
};

//-----------------------------------------------------------------------------
// class that manages file section parsers
class FileSectionMap : public map<string, FileSection*>
{
public:
	~FileSectionMap();
};

//-----------------------------------------------------------------------------
// Module Section
class FEBioModuleSection : public FileSection
{
public:
	FEBioModuleSection(FEFEBioImport* pim) : FileSection(pim) {}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Control Section
class FEBioControlSection : public FileSection
{
public:
	FEBioControlSection(FEFEBioImport* pim) : FileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	FESolver* BuildSolver(int nmod, FEM& fem);

	bool ParseCommonParams(XMLTag& tag);
	bool ParsePoroParams  (XMLTag& tag);
	bool ParseSoluteParams(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Material Section
class FEBioMaterialSection : public FileSection
{
public:
	FEBioMaterialSection(FEFEBioImport* pim) : FileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseMaterial					(XMLTag& tag, FEMaterial* pm);
	bool ParseElasticMaterial			(XMLTag& tag, FEElasticMaterial* pm);
	bool ParseTransIsoMaterial			(XMLTag& tag, FETransverselyIsotropic* pm);
	bool ParseRigidMaterial				(XMLTag& tag, FERigidMaterial* pm);
	bool ParseElasticMixture			(XMLTag& tag, FEElasticMixture* pm);
	bool ParseUncoupledElasticMixture	(XMLTag& tag, FEUncoupledElasticMixture* pm);
	bool ParseBiphasicMaterial			(XMLTag& tag, FEBiphasic* pm);
	bool ParseBiphasicSoluteMaterial	(XMLTag& tag, FEBiphasicSolute* pm);

protected:
	int	m_nmat;
};

//-----------------------------------------------------------------------------
// Geometry Section
class FEBioGeometrySection : public FileSection
{
private:
	enum {
		ET_HEX,
		ET_PENTA,
		ET_TET,
		ET_QUAD,
		ET_TRI,
		ET_TRUSS
	};
	
public:
	FEBioGeometrySection(FEFEBioImport* pim) : FileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseNodeSection       (XMLTag& tag);
	void ParseElementSection    (XMLTag& tag);
	void ParseElementDataSection(XMLTag& tag);
	void ParseGroupSection      (XMLTag& tag);

	void ReadSolidElement(XMLTag& tag, FESolidElement& el, int ntype, int nid, int gid, int nmat);
	void ReadShellElement(XMLTag& tag, FEShellElement& el, int ntype, int nid, int gid, int nmat);
	void ReadTrussElement(XMLTag& tag, FETrussElement& el, int ntype, int nid, int gid, int nmat);

	int ElementType(XMLTag& tag);
	int DomainType(int etype, FEMaterial* pmat);
	FEDomain* CreateDomain(int ntype, FEMesh* pm, FEMaterial* pmat);
};

//-----------------------------------------------------------------------------
// Boundary Section
class FEBioBoundarySection : public FileSection
{
public:
	FEBioBoundarySection(FEFEBioImport* pim) : FileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseBCFix               (XMLTag& tag);
	void ParseBCPrescribe         (XMLTag& tag);
	void ParseBCForce             (XMLTag& tag);
	void ParseBCPressure          (XMLTag& tag);
	void ParseBCTraction          (XMLTag& tag);
	void ParseBCPoroNormalTraction(XMLTag& tag);
	void ParseBCFluidFlux         (XMLTag& tag);
	void ParseBCSoluteFlux        (XMLTag &tag);
	void ParseBCHeatFlux          (XMLTag& tag);
	void ParseContactSection      (XMLTag& tag);
	void ParseConstraints         (XMLTag& tag);
	void ParseSpringSection       (XMLTag& tag);
	bool ParseSurfaceSection      (XMLTag& tag, FESurface& s, int nfmt);
};

//-----------------------------------------------------------------------------
// Initial Section
class FEBioInitialSection : public FileSection
{
public:
	FEBioInitialSection(FEFEBioImport* pim) : FileSection(pim){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Globals Section
class FEBioGlobalsSection : public FileSection
{
public:
	FEBioGlobalsSection(FEFEBioImport* pim) : FileSection(pim){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// LoadData Section
class FEBioLoadSection : public FileSection
{
public:
	FEBioLoadSection(FEFEBioImport* pim) : FileSection(pim) {}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Output Section
class FEBioOutputSection : public FileSection
{
public:
	FEBioOutputSection(FEFEBioImport* pim) : FileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseLogfile (XMLTag& tag);
	void ParsePlotfile(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Constraints Section
class FEBioConstraintsSection : public FileSection
{
public:
	FEBioConstraintsSection(FEFEBioImport* pim) : FileSection(pim){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Step Section
class FEBioStepSection : public FileSection
{
public:
	FEBioStepSection(FEFEBioImport* pim) : FileSection(pim){}
	void Parse(XMLTag& tag);
};

//=============================================================================
//! Implements a class to import FEBio input files
//!
class FEFEBioImport : public FEFileImport
{
public:
	// Element types
	enum { ET_HEX8, ET_PENTA6, ET_TET4, ET_UT4, ET_TETG1, ET_QUAD4, ET_TRI3, ET_TRUSS2 };

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

public:
	bool Load(FEM& fem, const char* szfile);

	FEM* GetFEM() { return m_pfem; }
	FEAnalysis*	GetStep() { return m_pStep; }

	int Version() { return m_nversion; }

protected:
	void ParseVersion			(XMLTag& tag);

public:
	FEM*		m_pfem;		//!< pointer to the fem class
	FEAnalysis*	m_pStep;	//!< pointer to current analysis step

protected:
	XMLReader	m_xml;	//!< the actual reader

public:
	int	m_ntet4;	// tetrahedral integration rule
	int m_nsteps;	// nr of step sections read
	int	m_nmat;		// nr of materials

protected:
	int	m_nversion;	// version of file
};
