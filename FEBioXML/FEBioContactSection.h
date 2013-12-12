#pragma once
#include "FEBioImport.h"
#include "FECore/FESurfacePairInteraction.h"

//-----------------------------------------------------------------------------
// Contact section (new in version 2.0)
class FEBioContactSection : public FEBioFileSection
{
public:
	FEBioContactSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseRigidWall            (XMLTag& tag);
	void ParseRigidInterface       (XMLTag& tag);
	void ParseRigidJoint           (XMLTag& tag);
	void ParseLinearConstraint     (XMLTag& tag);

protected:
	void ParseContactInterface(XMLTag& tag, FESurfacePairInteraction* pci);
	bool ParseSurfaceSection  (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);
};
