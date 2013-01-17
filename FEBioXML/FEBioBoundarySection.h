#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Boundary Section
class FEBioBoundarySection : public FEBioFileSection
{
public:
	FEBioBoundarySection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseBCFix               (XMLTag& tag);
	void ParseBCPrescribe         (XMLTag& tag);
	void ParseBCForce             (XMLTag& tag);
	void ParseBCPressure          (XMLTag& tag);
	void ParseBCTraction          (XMLTag& tag);
	void ParseBCPoroNormalTraction(XMLTag& tag);
	void ParseBCFluidFlux         (XMLTag& tag);
	void ParseBCSoluteFlux        (XMLTag &tag);
	void ParseBCHeatFlux          (XMLTag& tag);
	void ParseBCConvectiveHeatFlux(XMLTag& tag);
	void ParseContactSection      (XMLTag& tag);
	void ParseConstraints         (XMLTag& tag);
	void ParseSpringSection       (XMLTag& tag);
	bool ParseSurfaceSection      (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);

protected:
	void ParseRigidJoint      (XMLTag& tag);
	void ParseLinearConstraint(XMLTag& tag);
	void ParseRigidWall       (XMLTag& tag);
	void ParseRigidContact    (XMLTag& tag);
};
