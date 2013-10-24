#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Geometry Section
class FEBioGeometrySection : public FEBioFileSection
{
private:
	enum {
		ET_HEX8,
		ET_HEX20,
		ET_PENTA6,
		ET_TET4,
		ET_TET10,
		ET_QUAD4,
		ET_TRI3,
		ET_TRUSS2
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

	void ParseMesh(XMLTag& tag);

	void ReadSolidElement(XMLTag& tag, FESolidElement& el, int ntype, int nid, int nmat);
	void ReadShellElement(XMLTag& tag, FEShellElement& el, int ntype, int nid, int nmat);
	void ReadTrussElement(XMLTag& tag, FETrussElement& el, int ntype, int nid, int nmat);

	int ElementType(XMLTag& tag);
	int DomainType(int etype, FEMaterial* pmat);
	FEDomain* CreateDomain(int ntype, FEMesh* pm, FEMaterial* pmat);
};
