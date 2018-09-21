#pragma once
#include "FEBioImport.h"
#include "FEBModel.h"

//-----------------------------------------------------------------------------
// Geometry Section (base class)
class FEBioGeometrySection : public FEBioFileSection
{
public:
	FEBioGeometrySection(FEBioImport* pim) : FEBioFileSection(pim) {}

protected:
	void ReadElement(XMLTag& tag, FEElement& el, int nid);
};

//-----------------------------------------------------------------------------
class FEBioGeometrySection1x : public FEBioGeometrySection
{
protected:
	struct FEDOMAIN
	{
		FE_Element_Spec		elem;	// element type
		int					mat;	// material ID
		int					nel;	// number of elements
	};

public:
	FEBioGeometrySection1x(FEBioImport* pim) : FEBioGeometrySection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodeSection(XMLTag& tag);
	void ParseNodeSetSection(XMLTag& tag);
	void ParseElementSection(XMLTag& tag);
	void ParseElementDataSection(XMLTag& tag);
	void ParseElementData(FEElement& el, XMLTag& tag);
};

//-----------------------------------------------------------------------------
class FEBioGeometrySection2 : public FEBioGeometrySection
{
public:
	FEBioGeometrySection2(FEBioImport* pim) : FEBioGeometrySection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodeSection(XMLTag& tag);
	void ParseEdgeSection(XMLTag& tag);
	void ParseSurfaceSection   (XMLTag& tag);
	void ParseElementSection   (XMLTag& tag);
	void ParseNodeSetSection   (XMLTag& tag);
	void ParseElementSetSection(XMLTag& tag);
	void ParseElementDataSection(XMLTag& tag);
	void ParseElementData(FEElement& el, XMLTag& tag);
};

//-----------------------------------------------------------------------------
class FEBioGeometrySection25 : public FEBioGeometrySection
{
public:
	FEBioGeometrySection25(FEBioImport* pim) : FEBioGeometrySection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodeSection       (XMLTag& tag);
	void ParseDiscreteSetSection(XMLTag& tag);
	void ParseSurfacePairSection(XMLTag& tag);
	void ParseNodeSetPairSection(XMLTag& tag);
	void ParseNodeSetSetSection (XMLTag& tag);
	void ParsePartSection       (XMLTag& tag);
	void ParseInstanceSection   (XMLTag& tag);
	void ParseSurfaceSection    (XMLTag& tag);
	void ParseElementSection    (XMLTag& tag);
	void ParseNodeSetSection    (XMLTag& tag);
	void ParseEdgeSection       (XMLTag& tag);
	void ParseElementSetSection (XMLTag& tag);

	// New functions for parsing parts
	void ParsePart(XMLTag& tag, FEBModel::Part* part);
	void ParsePartNodeSection(XMLTag& tag, FEBModel::Part* part);
	void ParsePartElementSection(XMLTag& tag, FEBModel::Part* part);
	void ParsePartNodeSetSection(XMLTag& tag, FEBModel::Part* part);
	void ParsePartSurfaceSection(XMLTag& tag, FEBModel::Part* part);

protected:
	FEBModel			m_feb;
};

//-----------------------------------------------------------------------------
class FEBioGeometrySection3 : public FEBioGeometrySection
{
public:
	FEBioGeometrySection3(FEBioImport* pim) : FEBioGeometrySection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodeSection       (XMLTag& tag);
	void ParseDiscreteSetSection(XMLTag& tag);
	void ParseSurfacePairSection(XMLTag& tag);
	void ParseNodeSetPairSection(XMLTag& tag);
	void ParseNodeSetSetSection (XMLTag& tag);
	void ParsePartSection       (XMLTag& tag);
	void ParseInstanceSection   (XMLTag& tag);
	void ParseSurfaceSection    (XMLTag& tag);
	void ParseElementSection    (XMLTag& tag);
	void ParseNodeSetSection    (XMLTag& tag);
	void ParseEdgeSection       (XMLTag& tag);
	void ParseElementSetSection (XMLTag& tag);

	// New functions for parsing parts
	void ParsePart(XMLTag& tag, FEBModel::Part* part);
	void ParsePartNodeSection(XMLTag& tag, FEBModel::Part* part);
	void ParsePartElementSection(XMLTag& tag, FEBModel::Part* part);
	void ParsePartNodeSetSection(XMLTag& tag, FEBModel::Part* part);
	void ParsePartSurfaceSection(XMLTag& tag, FEBModel::Part* part);

protected:
	FEBModel			m_feb;
};
