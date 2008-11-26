#pragma once

#include "FileImport.h"
#include "XMLReader.h"

//-----------------------------------------------------------------------------
//! Implements a class to import FEBio input files

//! \todo: Define classes for each section.

class FEFEBioImport : public FEFileImport
{
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
	bool ParseControlSection    (XMLTag& tag);
	bool ParseMaterialSection   (XMLTag& tag);
	bool ParseGeometrySection   (XMLTag& tag);
	bool ParseNodeSection       (XMLTag& tag);
	bool ParseElementSection    (XMLTag& tag);
	bool ParseElementDataSection(XMLTag& tag);
	bool ParseGroupSection      (XMLTag& tag);
	bool ParseBoundarySection   (XMLTag& tag);
	bool ParseConstraints       (XMLTag& tag);
	bool ParseContactSection    (XMLTag& tag);
	bool ParseGlobalsSection    (XMLTag& tag);
	bool ParseLoadSection       (XMLTag& tag);
	bool ParseOutputSection     (XMLTag& tag);
	bool ParseStepSection       (XMLTag& tag);

	FEM*		m_pfem;		//!< pointer to the fem class
	FEAnalysis*	m_pStep;	//!< pointer to current analysis step

protected:
	XMLReader	m_xml;	//!< the actual reader

	int	m_nhex8;	// element type for hex8
	int m_nsteps;	// nr of step sections read
};
