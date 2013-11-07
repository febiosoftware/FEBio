#pragma once
#include "FEBioImport.h"
#include "FEBioMech/FETransverselyIsotropic.h"
#include "FEBioMech/FERigid.h"
#include "FEBioMech/FEElasticMixture.h"
#include "FEBioMech/FEUncoupledElasticMixture.h"
#include "FEBioMix/FEBiphasic.h"
#include "FEBioMix/FEBiphasicSolute.h"
#include "FEBioMix/FETriphasic.h"
#include "FEBioMix/FEMultiphasic.h"
#include "FEBioMech/FEViscoElasticMaterial.h"
#include "FEBioMech/FEUncoupledViscoElasticMaterial.h"
#include "FEBioMech/FEElasticMultigeneration.h"
#include "FEBioMech/FERemodelingElasticMaterial.h"

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
	bool ParseTransIsoMaterial			   (XMLTag& tag, FETransverselyIsotropic* pm);
	bool ParseReactionMaterial			   (XMLTag& tag, FEChemicalReaction* pm);
	bool ParseMultiphasicMaterial          (XMLTag& tag, FEMultiphasic* pm);
	bool ParseElasticMultigeneration	   (XMLTag &tag, FEElasticMultigeneration *pm);

	FEMaterial* CreateMaterial(XMLTag& tag);

protected:
	int	m_nmat;
};
