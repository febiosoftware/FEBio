#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Contact section (new in version 2.0)
class FEBioContactSection : public FEBioFileSection
{
public:
	FEBioContactSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseSlidingInterface     (XMLTag& tag);
	void ParseFacetSlidingInterface(XMLTag& tag);
	void ParseSlidingInterface2    (XMLTag& tag);
	void ParseSlidingInterface3    (XMLTag& tag);
	void ParseTiedInterface        (XMLTag& tag);
	void ParsePeriodicBoundary     (XMLTag& tag);
	void ParseSurfaceConstraint    (XMLTag& tag);
	void ParseRigidWall            (XMLTag& tag);
	void ParseRigidInterface       (XMLTag& tag);
	void ParseRigidJoint           (XMLTag& tag);
	void ParseLinearConstraint     (XMLTag& tag);

	bool ParseSurfaceSection      (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);
};
