#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Version 1.2
class FEBioLoadsSection1x : public FEFileSection
{
public:
	FEBioLoadsSection1x(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseNodalLoad  (XMLTag& tag);
	void ParseBodyForce  (XMLTag& tag);
	void ParseBodyLoad   (XMLTag& tag);
	void ParseSurfaceLoad(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Version 2.0
class FEBioLoadsSection2 : public FEFileSection
{
public:
	FEBioLoadsSection2(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseNodalLoad  (XMLTag& tag);
	void ParseBodyLoad   (XMLTag& tag);
	void ParseEdgeLoad   (XMLTag& tag);
	void ParseSurfaceLoad(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Version 2.5
class FEBioLoadsSection25 : public FEFileSection
{
public:
	FEBioLoadsSection25(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseNodalLoad  (XMLTag& tag);
	void ParseEdgeLoad   (XMLTag& tag);
	void ParseSurfaceLoad(XMLTag& tag);
	void ParseBodyLoad   (XMLTag& tag);
};
