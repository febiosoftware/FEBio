#pragma once
#include "FECoreBase.h"

class FEMaterialPoint;
class FEDomainMap;

class FEVectorGenerator : public FECoreBase
{
	DECLARE_SUPER_CLASS(FEVECTORGENERATOR_ID);

public:
	FEVectorGenerator(FEModel* fem);
	virtual ~FEVectorGenerator();

	// Call this function to calculate a vector at a material point
	virtual vec3d GetVector(const FEMaterialPoint& mp) = 0;
};

class FELocalVectorGenerator : public FEVectorGenerator
{
public:
	FELocalVectorGenerator(FEModel* fem);

	bool Init() override;

	vec3d GetVector(const FEMaterialPoint& mp) override;

protected:
	int	m_n[2];

	DECLARE_FECORE_CLASS();
};

class FEConstVectorGenerator : public FEVectorGenerator
{
public:
	FEConstVectorGenerator(FEModel* fem);

	bool Init() override;

	vec3d GetVector(const FEMaterialPoint& mp) override;

protected:
	vec3d	m_vec;

	DECLARE_FECORE_CLASS();
};

class FEUserVectorGenerator : public FEVectorGenerator
{
public:
	FEUserVectorGenerator(FEModel* fem);

	bool Init() override;

	vec3d GetVector(const FEMaterialPoint& mp) override;

	void SetData(FEDomainMap* map);

protected:
	FEDomainMap*	m_data;
};
