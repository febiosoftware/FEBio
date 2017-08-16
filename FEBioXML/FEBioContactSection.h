#pragma once
#include "FileImport.h"
#include <FECore/FESurfacePairInteraction.h>

//-----------------------------------------------------------------------------
// Contact section (new in version 2.0)
class FEBioContactSection : public FEFileSection
{
protected:
	//! missing slave surface
	class MissingSlaveSurface : public FEFileException
	{
	public:
		MissingSlaveSurface();
	};

	//! missing master surface
	class MissingMasterSurface : public FEFileException
	{
	public:
		MissingMasterSurface();
	};

public:
	FEBioContactSection(FEFileImport* pim) : FEFileSection(pim){}

protected:
	void ParseLinearConstraint     (XMLTag& tag);

protected:
	bool ParseSurfaceSection  (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);
};

//-----------------------------------------------------------------------------
// Version 2.0
class FEBioContactSection2 : public FEBioContactSection
{
public:
	FEBioContactSection2(FEFileImport* im) : FEBioContactSection(im){}
	void Parse(XMLTag& tag);

protected:
	void ParseRigidInterface(XMLTag& tag);
	void ParseRigidWall(XMLTag& tag);
	void ParseContactInterface(XMLTag& tag, FESurfacePairInteraction* pci);
};

//-----------------------------------------------------------------------------
// Version 2.5
class FEBioContactSection25 : public FEBioContactSection
{
public:
	FEBioContactSection25(FEFileImport* im) : FEBioContactSection(im){}
	void Parse(XMLTag& tag);

protected:
	void ParseRigidWall(XMLTag& tag);
	void ParseRigidSliding(XMLTag& tag);
	void ParseContactInterface(XMLTag& tag, FESurfacePairInteraction* pci);
};
