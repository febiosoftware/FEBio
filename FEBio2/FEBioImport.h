#pragma once

#include "FileImport.h"
#include "FECore/XMLReader.h"
#include "FEBioLib/FETransverselyIsotropic.h"
#include "FEBioLib/FERigid.h"
#include "FEBioLib/FEElasticMixture.h"
#include "FEBioLib/FEUncoupledElasticMixture.h"
#include "FEBioLib/FEBiphasic.h"
#include "FEBioLib/FEBiphasicSolute.h"
#include "FEBioLib/FETriphasic.h"
#include "FEAnalysisStep.h"
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

	FEM* GetFEM();
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

//-----------------------------------------------------------------------------
// Import section
class FEBioImportSection : public FEBioFileSection
{
public:
	FEBioImportSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Module Section
class FEBioModuleSection : public FEBioFileSection
{
public:
	FEBioModuleSection(FEFEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Control Section
class FEBioControlSection : public FEBioFileSection
{
public:
	FEBioControlSection(FEFEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	FESolver* BuildSolver(int nmod, FEM& fem);

	bool ParseCommonParams	   (XMLTag& tag);
	void ParseSolidParams	   (XMLTag& tag);
	void ParsePoroParams	   (XMLTag& tag);
	void ParseSoluteParams	   (XMLTag& tag);
	void ParseTriphasicParams  (XMLTag& tag);
	void ParseLinearSolidParams(XMLTag& tag);
	void ParseHeatParams       (XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Material Section
class FEBioMaterialSection : public FEBioFileSection
{
public:
	FEBioMaterialSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
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
	bool ParseSoluteMaterial			(XMLTag& tag, FESolute* pm);
	bool ParseTriphasicMaterial			(XMLTag& tag, FETriphasic* pm);
	bool ParseNestedMaterial			(XMLTag& tag, FENestedMaterial* pm);

protected:
	int	m_nmat;
};

//-----------------------------------------------------------------------------
// Geometry Section
class FEBioGeometrySection : public FEBioFileSection
{
private:
	enum {
		ET_HEX,
		ET_HEX20,
		ET_PENTA,
		ET_TET,
		ET_QUAD,
		ET_TRI,
		ET_TRUSS
	};

	struct FEDOMAIN 
	{
		int		mat;	// material ID
		int		elem;	// element type
		int		nel;	// number of elements
	};
	
public:
	FEBioGeometrySection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseNodeSection       (XMLTag& tag);
	void ParseElementSection    (XMLTag& tag);
	void ParseElementDataSection(XMLTag& tag);
	void ParseNodeSetSection    (XMLTag& tag);
	void ParsePartSection       (XMLTag& tag);

	void ReadSolidElement(XMLTag& tag, FESolidElement& el, int ntype, int nid, int nmat);
	void ReadShellElement(XMLTag& tag, FEShellElement& el, int ntype, int nid, int nmat);
	void ReadTrussElement(XMLTag& tag, FETrussElement& el, int ntype, int nid, int nmat);

	int ElementType(XMLTag& tag);
	int DomainType(int etype, FEMaterial* pmat);
	FEDomain* CreateDomain(int ntype, FEMesh* pm, FEMaterial* pmat);
};

//-----------------------------------------------------------------------------
// Boundary Section
class FEBioBoundarySection : public FEBioFileSection
{
public:
	FEBioBoundarySection(FEFEBioImport* pim) : FEBioFileSection(pim){}
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
// Loads Section (new in version 1.2)
class FEBioLoadsSection : public FEBioBoundarySection
{
public:
	FEBioLoadsSection(FEFEBioImport* pim) : FEBioBoundarySection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseBodyForce(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Contact section (new in version 2.0)
class FEBioContactSection : public FEBioFileSection
{
public:
	FEBioContactSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseSlidingInterface     (XMLTag& tag);
	void ParseFacetSlidingInterface(XMLTag& tag);
	void ParseSlidingInterface2    (XMLTag& tag);
	void ParseSlidingInterface3    (XMLTag& tag);
	void ParseTiedInterface        (XMLTag& tag);
	void ParsePeriodicBoundary     (XMLTag& tag);
	void ParseSurfaceConstraint    (XMLTag& tag);
	void ParseRigidWall            (XMLTag& tag);
	void ParseRigidInterface       (XMLTag& tag);
	void ParseRigidJoint           (XMLTag& tag);
	void ParseLinearConstraint     (XMLTag& tag);

	bool ParseSurfaceSection      (XMLTag& tag, FESurface& s, int nfmt);
};

//-----------------------------------------------------------------------------
// Initial Section
class FEBioInitialSection : public FEBioFileSection
{
public:
	FEBioInitialSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Globals Section
class FEBioGlobalsSection : public FEBioFileSection
{
public:
	FEBioGlobalsSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse            (XMLTag& tag);

protected:
	void ParseBodyForce   (XMLTag& tag);	// only for versions < 1.2
	void ParseConstants   (XMLTag& tag);
	void ParseGSSoluteData(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// LoadData Section
class FEBioLoadDataSection : public FEBioFileSection
{
public:
	FEBioLoadDataSection(FEFEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Output Section
class FEBioOutputSection : public FEBioFileSection
{
public:
	FEBioOutputSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseLogfile (XMLTag& tag);
	void ParsePlotfile(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Constraints Section
class FEBioConstraintsSection : public FEBioFileSection
{
public:
	FEBioConstraintsSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseRigidConstraint(XMLTag& tag);
	void ParsePointConstraint(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Step Section
class FEBioStepSection : public FEBioFileSection
{
public:
	FEBioStepSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);
};

//=============================================================================
//! Implements a class to import FEBio input files
//!
class FEFEBioImport : public FEFileImport
{
public:
	// Element types
	enum { ET_HEX8, ET_HEX20, ET_PENTA6, ET_TET4, ET_UT4, ET_TETG1, ET_QUAD4, ET_TRI3, ET_TRUSS2 };

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
	bool Load(FEM& fem, const char* szfile);

	FEM* GetFEM() { return m_pfem; }
	FEAnalysis*	GetStep() { return m_pStep; }

	int Version() { return m_nversion; }

	bool ReadParameter(XMLTag& tag, FEParameterList& pl);

	void ReadList(XMLTag& tag, vector<int>& l);

protected:
	void ParseVersion			(XMLTag& tag);

public:
	FEM*		m_pfem;		//!< pointer to the fem class
	FEAnalysis*	m_pStep;	//!< pointer to current analysis step

public:
	int	m_ntet4;	// tetrahedral integration rule
	int	m_nut4;		// integration rule for stabilization of UT4
	int m_nsteps;	// nr of step sections read
	int	m_nmat;		// nr of materials

	bool	m_b3field;	// three-field element flag
	int		m_nhex8;	// hex integration rule

protected:
	int	m_nversion;	// version of file
};
