#pragma once
#include "FEBioImport.h"
#include "FECore/FESurfacePairInteraction.h"

//-----------------------------------------------------------------------------
// Boundary Section
class FEBioBoundarySection : public FEFileSection
{
public:
	FEBioBoundarySection(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseBCFix         (XMLTag& tag);
	void ParseBCFix20       (XMLTag& tag);
	void ParseBCFix25       (XMLTag& tag);
	void ParseBCPrescribe   (XMLTag& tag);
	void ParseBCPrescribe20 (XMLTag& tag);
	void ParseBCPrescribe25 (XMLTag& tag);
	void ParseContactSection(XMLTag& tag);
	void ParseConstraints   (XMLTag& tag);
	void ParseSpringSection (XMLTag& tag);
	void ParseBCRigid       (XMLTag& tag);	// version 2.5 and up
	void ParseRigidBody     (XMLTag& tag);  // version 2.5 and up
	void ParseBC            (XMLTag& tag);  // version 2.5 and up

	void ParsePeriodicLinearConstraint  (XMLTag& tag); // version 2.5 (temporary construction)
	void ParsePeriodicLinearConstraint2O(XMLTag& tag); // version 2.5 (temporary construction)
	void ParseMergeConstraint(XMLTag& tag); // version 2.5

protected:
	void ParseContactInterface(XMLTag& tag, FESurfacePairInteraction* psi);
	bool ParseSurfaceSection  (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);

protected:
	void ParseRigidJoint      (XMLTag& tag);
	void ParseLinearConstraint(XMLTag& tag);
	void ParseRigidWall       (XMLTag& tag);
	void ParseRigidContact    (XMLTag& tag);
};
