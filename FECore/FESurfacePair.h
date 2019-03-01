#pragma once
#include "fecore_api.h"
#include <string>

//-----------------------------------------------------------------------------
class FEMesh;
class FEFacetSet;
class DumpStream;

//-----------------------------------------------------------------------------
class FECORE_API FESurfacePair
{
public:
	FESurfacePair(FEMesh* pm);

	void SetName(const std::string& name);
	const std::string& GetName() const;

	FEFacetSet* GetMasterSurface();
	void SetMasterSurface(FEFacetSet* pf);

	FEFacetSet* GetSlaveSurface();
	void SetSlaveSurface(FEFacetSet* pf);

	void Serialize(DumpStream& ar);

private:
	std::string		m_name;
	FEFacetSet*		m_master;
	FEFacetSet*		m_slave;
	FEMesh*			m_mesh;
};
