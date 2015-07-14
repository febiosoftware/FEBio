#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Geometry Section
class FEBioGeometrySection : public FEBioFileSection
{
private:
	struct FEDOMAIN 
	{
		FE_Element_Shape	elem;	// element type
		int					mat;	// material ID
		int					nel;	// number of elements
	};
	
public:
	FEBioGeometrySection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseNodeSection       (XMLTag& tag);
	void ParseElementSection    (XMLTag& tag);
	void ParseElementSection20  (XMLTag& tag);
	void ParseElementDataSection(XMLTag& tag);
	void ParseNodeSetSection    (XMLTag& tag);
	void ParseSurfaceSection    (XMLTag& tag);
	void ParseElementSetSection (XMLTag& tag);

	void ParseMesh(XMLTag& tag);

	void ReadSolidElement(XMLTag& tag, FESolidElement& el, int ntype, int nid, int nmat);
	void ReadShellElement(XMLTag& tag, FEShellElement& el, int ntype, int nid, int nmat);
	void ReadTrussElement(XMLTag& tag, FETrussElement& el, int ntype, int nid, int nmat);

	FE_Element_Shape ElementShape(XMLTag& tag);
	int DomainType(FE_Element_Shape eshape, FEMaterial* pmat);
	FEDomain* CreateDomain(int ntype, FEMesh* pm, FEMaterial* pmat);


protected:
	vector<FEDOMAIN>	m_dom;
};
