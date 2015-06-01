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
	void ParseRigidSphericalJoint  (XMLTag& tag);
    void ParseRigidRevoluteJoint   (XMLTag& tag);
    void ParseRigidPrismaticJoint  (XMLTag& tag);
    void ParseRigidCylindricalJoint(XMLTag& tag);
    void ParseRigidPlanarJoint     (XMLTag& tag);
    void ParseRigidSpring          (XMLTag& tag);
    void ParseRigidDamper          (XMLTag& tag);
    void ParseRigidAngularDamper   (XMLTag& tag);
    void ParseRigidContractileForce(XMLTag& tag);
	void ParseLinearConstraint     (XMLTag& tag);

protected:
	void ParseContactInterface(XMLTag& tag, FESurfacePairInteraction* pci);
	bool ParseSurfaceSection  (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);

protected:
	bool BuildSurface(FESurface& s, FEFacetSet& f, bool bnodal);
};
