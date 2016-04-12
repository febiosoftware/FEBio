#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
class FEFaceSet;
class FEEdge;
class FEFacetSet;
class FESegmentSet;

//-----------------------------------------------------------------------------
// Loads Section (new in version 1.2)
class FEBioLoadsSection : public FEBioFileSection
{
public:
	FEBioLoadsSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseNodalLoad    (XMLTag& tag);
	void ParseNodalLoad25  (XMLTag& tag);
	void ParseBodyForce    (XMLTag& tag);
	void ParseBodyLoad     (XMLTag& tag);
	void ParseBodyLoad20   (XMLTag& tag);
	void ParseEdgeLoad     (XMLTag& tag);
	void ParseSurfaceLoad  (XMLTag& tag);
	void ParseSurfaceLoad20(XMLTag& tag);
	void ParseSurfaceLoad25(XMLTag& tag);

protected:
	bool BuildEdge   (FEEdge&    s, FESegmentSet& f);
};
