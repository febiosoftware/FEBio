#pragma once
#include "FEBioImport.h"
#include "FEBioMech/FEContactInterface.h"

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
	void ParseContactSection      (XMLTag& tag);
	void ParseConstraints         (XMLTag& tag);
	void ParseSpringSection       (XMLTag& tag);

protected:
	void ParseContactInterface(XMLTag& tag, FEContactInterface* pci);
	bool ParseSurfaceSection      (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);

protected:
	void ParseRigidJoint      (XMLTag& tag);
	void ParseLinearConstraint(XMLTag& tag);
	void ParseRigidWall       (XMLTag& tag);
	void ParseRigidContact    (XMLTag& tag);
};
