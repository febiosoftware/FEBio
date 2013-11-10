#pragma once
#include "FEBioImport.h"
#include "FEBioMech/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
// Material Section
class FEBioMaterialSection : public FEBioFileSection
{
public:
	FEBioMaterialSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseMaterial					   (XMLTag& tag, FEMaterial* pm);
	bool ParseElasticMaterial			   (XMLTag& tag, FEElasticMaterial* pm);

	FEMaterial* CreateMaterial(XMLTag& tag);

protected:
	int	m_nmat;
};
