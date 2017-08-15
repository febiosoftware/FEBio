#pragma once
#include "FEBioImport.h"
#include "FEBModel.h"

//-----------------------------------------------------------------------------
// Geometry Section (base class)
class FEBioGeometrySection : public FEBioFileSection
{
private:
	struct FEDOMAIN 
	{
		FE_Element_Spec		elem;	// element type
		int					mat;	// material ID
		int					nel;	// number of elements
	};

public:
	FEBioGeometrySection(FEBioImport* pim) : FEBioFileSection(pim) {}

protected:
	void ParseNodeSection       (XMLTag& tag);
	void ParseElementSection(XMLTag& tag);
	void ParseElementSection20  (XMLTag& tag);
	void ParseElementDataSection(XMLTag& tag);
	void ParseNodeSetSection    (XMLTag& tag);
	void ParseDiscreteSetSection(XMLTag& tag);
	void ParseEdgeSection       (XMLTag& tag);
	void ParseSurfaceSection    (XMLTag& tag);
	void ParseNodeSetPairSection(XMLTag& tag);
	void ParseNodeSetSetSection (XMLTag& tag);
	void ParseSurfacePairSection(XMLTag& tag);
	void ParseElementSetSection (XMLTag& tag);
	void ParsePartSection       (XMLTag& tag);
	void ParseInstanceSection   (XMLTag& tag);
	void ParseElementData(FEElement& el, XMLTag& tag);

	// New functions for parsing parts
	void ParsePart(XMLTag& tag, FEBModel::Part* part);
	void ParseNodeSection25   (XMLTag& tag, FEBModel::Part* part);
	void ParseElementSection25(XMLTag& tag, FEBModel::Part* part);
	void ParseNodeSetSection25(XMLTag& tag, FEBModel::Part* part);
	void ParseSurfaceSection25(XMLTag& tag, FEBModel::Part* part);

	void ParseMesh(XMLTag& tag);

	void ReadElement(XMLTag& tag, FEElement& el, int nid);

	FE_Element_Spec ElementSpec(const char* sz);
	FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);

protected:
	vector<FEDOMAIN>	m_dom;
	FEBModel			m_feb;	// User by 2.5 when defining the geometry in parts
};

//-----------------------------------------------------------------------------
class FEBioGeometrySection1x : public FEBioGeometrySection
{
public:
	FEBioGeometrySection1x(FEBioImport* pim) : FEBioGeometrySection(pim) {}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
class FEBioGeometrySection2 : public FEBioGeometrySection
{
public:
	FEBioGeometrySection2(FEBioImport* pim) : FEBioGeometrySection(pim) {}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
class FEBioGeometrySection25 : public FEBioGeometrySection
{
public:
	FEBioGeometrySection25(FEBioImport* pim) : FEBioGeometrySection(pim) {}
	void Parse(XMLTag& tag);
};
