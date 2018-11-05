#pragma once
#include "FEDataGenerator.h"
#include "FEFunction1D.h"
#include "FEClosestPointProjection.h"
#include "FENormalProjection.h"

class FEModel;
class FESurface;

class FESurfaceToSurfaceMap : public FEDataGenerator
{
public:
	FESurfaceToSurfaceMap(FEModel* fem);
	~FESurfaceToSurfaceMap();

	bool Init() override;

	void value(const vec3d& x, double& data) override;

private:
	std::string		m_surfName1;
	std::string		m_surfName2;
	FEFunction1D*	m_func;
	
private:
	FESurface*	m_surf1;
	FESurface*	m_surf2;
	FEClosestPointProjection*	m_ccp;
	FENormalProjection*			m_npr;

	DECLARE_FECORE_CLASS();
};
