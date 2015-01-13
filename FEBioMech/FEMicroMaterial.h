#pragma once
#include "FEElasticMaterial.h"
#include "FECore/FEModel.h"
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Material point class for the micro-material
class FEMicroMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEMicroMaterialPoint(FEMaterialPoint* mp);

	//! Initialize material point data
	void Init(bool bflag);

	//! create a shallow copy
	FEMaterialPoint* Copy();

	//! serialize material point data
	void Serialize(DumpFile& ar);

	//! stream material point data
	void ShallowCopy(DumpStream& dmp, bool bsave);

public:
	tens4ds	m_Ka;	//!< averaged material stiffness
};

//-----------------------------------------------------------------------------
//! The micro-material implements material homogenization. The stress and tangents
//! are calculated by solving a micro-structural RVE problem and return the
//! averaged stress and condensed tangents.
//!
class FEMicroMaterial :	public FEElasticMaterial
{
public:
	FEMicroMaterial(FEModel* pfem);
	~FEMicroMaterial(void);

public:
	char	m_szrve[256];	//!< filename for RVE file

protected:
	FEModel		m_rve;	//!< the RVE (Representive Volume Element)
	bool		m_brve;	//!< flag indicating whether RVE was read in
	double		m_V0;	//!< initial volume of RVE

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization
	void Init();

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData();

protected:
	void PrepRVE();
	bool SolveRVE(mat3d& F);
	void UpdateBC(mat3d& F);
	mat3ds AveragedStress(FEMaterialPoint& pt);
	tens4ds AveragedStiffness(FEMaterialPoint& pt);

public:
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
