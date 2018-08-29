#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
class FEBioDiscreteSection : public FEFileSection
{
public:
	FEBioDiscreteSection(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseSpringSection  (XMLTag& tag);
	void ParseRigidAxialForce(XMLTag& tag);
};

//-----------------------------------------------------------------------------
class FEBioDiscreteSection25 : public FEFileSection
{
public:
	FEBioDiscreteSection25(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseRigidAxialForce(XMLTag& tag);
	void ParseRigidCable(XMLTag& tag);
};
