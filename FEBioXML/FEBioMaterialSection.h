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
	bool ParseElasticMixture			   (XMLTag& tag, FEElasticMixture* pm);
	bool ParseUncoupledElasticMixture	   (XMLTag& tag, FEUncoupledElasticMixture* pm);
	bool ParseBiphasicMaterial		  	   (XMLTag& tag, FEBiphasic* pm);
	bool ParseBiphasicSoluteMaterial  	   (XMLTag& tag, FEBiphasicSolute* pm);
	bool ParseSoluteMaterial			   (XMLTag& tag, FESolute* pm);
	bool ParseReactionMaterial			   (XMLTag& tag, FEChemicalReaction* pm);
	bool ParseTriphasicMaterial  		   (XMLTag& tag, FETriphasic* pm);
	bool ParseMultiphasicMaterial          (XMLTag& tag, FEMultiphasic* pm);
	bool ParseViscoElasticMaterial		   (XMLTag& tag, FEViscoElasticMaterial* pm);
	bool ParseUncoupledViscoElasticMaterial(XMLTag& tag, FEUncoupledViscoElasticMaterial* pm);
	bool ParseElasticMultigeneration	   (XMLTag &tag, FEElasticMultigeneration *pm);
	bool ParseRemodelingSolid			   (XMLTag &tag, FERemodelingElasticMaterial *pm);

protected:
	int	m_nmat;
};
