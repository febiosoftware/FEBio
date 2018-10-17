#pragma once
#include "FECoreBase.h"

class FEModel;
class FEMaterialPoint;
class FEDomainMap;

//-----------------------------------------------------------------------------
// Base class for vector generators, which generate a vec3d value at a matererial point
class FEVectorGenerator : public FECoreBase
{
	DECLARE_SUPER_CLASS(FEVECTORGENERATOR_ID);

public:
	FEVectorGenerator(FEModel* fem);
	virtual ~FEVectorGenerator();

	// Call this function to calculate a vector at a material point
	virtual vec3d GetVector(const FEMaterialPoint& mp) = 0;
};

//-----------------------------------------------------------------------------
// This class calculates a vector based on local element node numbering
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

//-----------------------------------------------------------------------------
// This generates a constant vector
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

//-----------------------------------------------------------------------------
class FECylindricalVectorGenerator: public FEVectorGenerator
{
public:
	FECylindricalVectorGenerator(FEModel* fem);

	bool Init() override;

	vec3d GetVector(const FEMaterialPoint& mp) override;

protected:
	vec3d	m_center;		// center of map
	vec3d	m_axis;			// cylinder axis
	vec3d	m_vector;		// reference direction

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FESphericalVectorGenerator : public FEVectorGenerator
{
public:
	FESphericalVectorGenerator(FEModel* fem);

	bool Init() override;

	vec3d GetVector(const FEMaterialPoint& mp) override;

protected:
	vec3d	m_center;
	vec3d	m_vector;

	DECLARE_FECORE_CLASS();
};


//-----------------------------------------------------------------------------
// This reads a vector from a domain map
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
