#pragma once
#include "FEElasticMaterial.h"
#include "FECore/FEModel.h"
#include "FECore/FEMaterial.h"
#include "FEBioMech/FEPeriodicBoundary1O.h"

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
	mat3d      m_PK1;			// LTE - 1st Piola-Kirchhoff stress
	mat3ds     m_S;				// LTE - 2nd Piola-Kirchhoff stress
	mat3d      m_PK1_prev;

	mat3ds     m_inf_str;		// LTE - infinitesimal strain
	mat3ds     m_E;				// LTE - Green-Lagrange strain
	mat3ds     m_e;				// LTE - Euler-Almansi strain
	
	double     m_macro_energy;	// LTE - Macroscopic strain energy
	double	   m_micro_energy;	// LTE - Volume-average of strain energy throughout the RVE solution
	double	   m_energy_diff;	// LTE - Difference between macro energy and volume averaged energy of RVE (should be zero) 

	double	   m_macro_energy_inc;	// LTE - Macroscopic strain energy increment
	double	   m_micro_energy_inc;	// LTE - Microscopic strain energy increment

	bool    m_rve_init;			// LTE - Flag indicating that the rve has been initialized
	FEModel m_rve;				// LTE - Current copy of the rve		
	FEModel m_rve_prev;			// LTE - Previous converged state of the rve

	tens4ds	m_Ka;				// LTE - Averaged rank 4 material stiffness
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
	char	m_szbc[256];	//!< name of nodeset defining boundary
	bool	m_bperiodic;	//!< periodic bc flag
	FEModel	m_mrve;			//!< the master RVE (Representive Volume Element)

protected:
	double	m_V0;			//!< initial volume of RVE
	vector<int> m_BN;		//!< boundary node flags

	double m_bb_x; double m_bb_y; double m_bb_z;  // LTE - RVE bounding box
	int m_num_ext_node;

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	virtual mat3ds Stress1O(FEMaterialPoint& pt, int plot_on, int int_pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization
	bool Init();

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData();

protected:
	bool PrepRVE();
	bool PrepDisplacementBC();
	bool PrepPeriodicBC();
	void FindBoundaryNodes();

	double InitVolume();

	void UpdateBC(FEModel& rve, mat3d& F);
	
	mat3ds AveragedStress(FEModel& rve, FEMaterialPoint& pt);
	tens4ds AveragedStiffness(FEModel& rve, FEMaterialPoint& pt);

	mat3d AveragedStressPK1(FEModel& rve, FEMaterialPoint &mp);
	mat3ds AveragedStressPK2(FEModel& rve, FEMaterialPoint &mp);

	void calc_energy_diff(FEMaterialPoint& pt);

public:
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
