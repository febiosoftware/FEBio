#pragma once

#include "FileImport.h"
#include "XMLReader.h"

//-----------------------------------------------------------------------------
//! Implements a class to import FEBio input files

//! \todo: Define classes for each section.

class FEFEBioImport : public FEFileImport
{
	// Element types
	enum { ET_HEX8, ET_PENTA6, ET_TET4, ET_QUAD4, ET_TRI3, ET_TRUSS2 };

	// element classes
	enum { EC_STRUCT, EC_RIGID, EC_PORO };

public:
	class InvalidVersion{};
	class InvalidMaterial
	{ 
	public: 
		InvalidMaterial(int nel) : m_nel(nel){}
		int m_nel; 
	};

public:
	bool Load(FEM& fem, const char* szfile);

protected:
	bool ParseModuleSection     (XMLTag& tag);
	bool ParseControlSection    (XMLTag& tag);
	bool ParseMaterialSection   (XMLTag& tag);
	bool ParseGeometrySection   (XMLTag& tag);
	bool ParseNodeSection       (XMLTag& tag);
	bool ParseElementSection    (XMLTag& tag);
	bool ParseElementDataSection(XMLTag& tag);
	bool ParseGroupSection      (XMLTag& tag);
	bool ParseBoundarySection   (XMLTag& tag);
	bool ParseInitialSection    (XMLTag& tag);
	bool ParseConstraints       (XMLTag& tag);
	bool ParseContactSection    (XMLTag& tag);
	bool ParseGlobalsSection    (XMLTag& tag);
	bool ParseLoadSection       (XMLTag& tag);
	bool ParseOutputSection     (XMLTag& tag);
	bool ParseStepSection       (XMLTag& tag);
	bool ParseSurfaceSection    (XMLTag& tag, FESurface& s, int nfmt);

	void ReadSolidElement(XMLTag& tag, FESolidElement& el, int ntype, int nid, int gid, int nmat);
	void ReadShellElement(XMLTag& tag, FEShellElement& el, int ntype, int nid, int gid, int nmat);
	void ReadTrussElement(XMLTag& tag, FETrussElement& el, int ntype, int nid, int gid, int nmat);

	FEM*		m_pfem;		//!< pointer to the fem class
	FEAnalysis*	m_pStep;	//!< pointer to current analysis step

protected:
	XMLReader	m_xml;	//!< the actual reader

	int m_nsteps;	// nr of step sections read
};
