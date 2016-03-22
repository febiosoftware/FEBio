#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// MeshData Section (introduced in febio_spec 2.5)
class FEBioMeshDataSection : public FEBioFileSection
{
	struct ELEMENT_DATA
	{
		int		nid;
		int		nval;	// number of values read
		double	val[FEElement::MAX_NODES];	// scalar value
	};

public:
	FEBioMeshDataSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseShellThickness(XMLTag& tag);
	void ParseMaterialFibers(XMLTag& tag);
	void ParseMaterialAxes  (XMLTag& tag);
	void ParseMaterialData  (XMLTag& tag, const string& name);

private:
	void ParseElementData(XMLTag& tag, vector<ELEMENT_DATA>& values, int nvalues);

private:
	vector<FEElement*> m_pelem;
};
