#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
class FEBioDiscreteSection : public FEBioFileSection
{
public:
	FEBioDiscreteSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseDiscreteSection(XMLTag& tag);
	void ParseDiscreteSection25(XMLTag& tag);

protected:
	void ParseSpringSection  (XMLTag& tag);
	void ParseRigidAxialForce(XMLTag& tag);
};
