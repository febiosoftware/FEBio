#pragma once
#include "FEBioImport.h"
#include "FECore/FESurfacePairConstraint.h"
#include <map>

//-----------------------------------------------------------------------------
// Boundary Section
class FEBioBoundarySection : public FEFileSection
{
public:
	FEBioBoundarySection(FEFileImport* pim) : FEFileSection(pim){}

protected:
	void ParseBCFix         (XMLTag& tag);
	void ParseBCPrescribe   (XMLTag& tag);
	void ParseContactSection(XMLTag& tag);
	void ParseConstraints   (XMLTag& tag);
	void ParseSpringSection (XMLTag& tag);

protected:
	void ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* psi);
	bool ParseSurfaceSection  (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);

protected:
	void ParseRigidJoint      (XMLTag& tag);
	void ParseLinearConstraint(XMLTag& tag);
	void ParseRigidWall       (XMLTag& tag);
	void ParseRigidContact    (XMLTag& tag);
};

//-----------------------------------------------------------------------------
// older formats (1.2)
class FEBioBoundarySection1x : public FEBioBoundarySection
{
public:
	FEBioBoundarySection1x(FEFileImport* imp) : FEBioBoundarySection(imp){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// version 2.0
class FEBioBoundarySection2 : public FEBioBoundarySection
{
public:
	FEBioBoundarySection2(FEFileImport* imp) : FEBioBoundarySection(imp){}
	void Parse(XMLTag& tag);

protected:
	void ParseBCFix      (XMLTag& tag);
	void ParseBCPrescribe(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// version 2.5
class FEBioBoundarySection25 : public FEBioBoundarySection
{
public:
	FEBioBoundarySection25(FEFileImport* imp) : FEBioBoundarySection(imp){}
	void Parse(XMLTag& tag);

protected:
	void ParseBCFix      (XMLTag& tag);
	void ParseBCPrescribe(XMLTag& tag);
	void ParseBCRigid    (XMLTag& tag);
	void ParseRigidBody  (XMLTag& tag);
	void ParseBC         (XMLTag& tag);

	void ParsePeriodicLinearConstraint  (XMLTag& tag); // version 2.5 (temporary construction)
	void ParsePeriodicLinearConstraint2O(XMLTag& tag); // version 2.5 (temporary construction)
	void ParseMergeConstraint           (XMLTag& tag); // version 2.5

protected:
	void BuildNodeSetMap();

private:
	std::map<std::string, FENodeSet*>	m_NodeSet;	// map for faster lookup of node sets
};
