#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Material Section
class FEBioMaterialSection : public FEBioFileSection
{
public:
	FEBioMaterialSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseMaterial  (XMLTag& tag, FEMaterial* pm);
	bool ParseFiberTag  (XMLTag& tag, FEMaterial* pm);
	bool ParseMatAxisTag(XMLTag& tag, FEMaterial* pm);

	FEMaterial* CreateMaterial(XMLTag& tag);

protected:
	int	m_nmat;
};
