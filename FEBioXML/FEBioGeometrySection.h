#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Geometry Section
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
	FEBioGeometrySection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void Parse12(XMLTag& tag);
	void Parse20(XMLTag& tag);
	void Parse25(XMLTag& tag);

protected:
	void ParseNodeSection       (XMLTag& tag, FEBioImport::Part* part = 0);
	void ParseElementSection    (XMLTag& tag);
	void ParseElementSection20  (XMLTag& tag, FEBioImport::Part* part = 0);
	void ParseElementDataSection(XMLTag& tag);
	void ParseNodeSetSection    (XMLTag& tag, FEBioImport::Part* part = 0);
	void ParseDiscreteSetSection(XMLTag& tag);
	void ParseEdgeSection       (XMLTag& tag);
	void ParseSurfaceSection    (XMLTag& tag, FEBioImport::Part* part = 0);
	void ParseNodeSetPairSection(XMLTag& tag);
	void ParseSurfacePairSection(XMLTag& tag);
	void ParseElementSetSection (XMLTag& tag);
	void ParsePartSection       (XMLTag& tag);
	void ParseElementData(FEElement& el, XMLTag& tag);

	void ParseMesh(XMLTag& tag);

	void ReadElement(XMLTag& tag, FEElement& el, int nid);

	void ParsePart(XMLTag& tag, FEBioImport::Part* part);

	FE_Element_Spec ElementSpec(const char* sz);
	FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);

protected:
	vector<FEDOMAIN>	m_dom;
};
