#pragma once
#include "FEBioImport.h"
#include "FECore/FESurfacePairInteraction.h"

//-----------------------------------------------------------------------------
// Boundary Section
class FEBioBoundarySection : public FEBioFileSection
{
public:
	FEBioBoundarySection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseBCFix               (XMLTag& tag);
	void ParseBCFix20             (XMLTag& tag);
	void ParseBCPrescribe         (XMLTag& tag);
	void ParseBCPrescribe20       (XMLTag& tag);
	void ParseContactSection      (XMLTag& tag);
	void ParseConstraints         (XMLTag& tag);
	void ParseSpringSection       (XMLTag& tag);

protected:
	void ParseContactInterface(XMLTag& tag, FESurfacePairInteraction* psi);
	bool ParseSurfaceSection  (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);

protected:
	void ParseRigidJoint      (XMLTag& tag);
	void ParseLinearConstraint(XMLTag& tag);
	void ParseRigidWall       (XMLTag& tag);
	void ParseRigidContact    (XMLTag& tag);
};
