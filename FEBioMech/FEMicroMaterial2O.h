#pragma once
#include "FEElasticMaterial.h"
#include "FECore/FEModel.h"
#include "FECore/FEMaterial.h"
#include "FEBioMech/FEPeriodicBoundary2O.h"
#include "FECore/tens3d.h"
#include "FECore/tens4d.h"
#include "FECore/tens5d.h"
#include "FECore/tens6d.h"

//-----------------------------------------------------------------------------
//! Material point class for the micro-material
class FEMicroMaterialPoint2O : public FEMaterialPoint
{
public:
	//! constructor
	FEMicroMaterialPoint2O(FEMaterialPoint* mp);

	//! Initialize material point data
	void Init(bool bflag);

	//! create a shallow copy
	FEMaterialPoint* Copy();

	//! serialize material point data
	void Serialize(DumpFile& ar);

	//! stream material point data
	void ShallowCopy(DumpStream& dmp, bool bsave);

public:
	tens3drs   m_G;				// LTE - Deformation Hessian
	tens3ds    m_tau;			// LTE - Cauchy stress moment
	
	mat3d      m_PK1;			// LTE - 1st Piola-Kirchhoff stress
	tens3drs   m_QK1;			// LTE - 1st Piola-Kirchhoff stress moment

	mat3ds     m_S;				// LTE - 2nd Piola-Kirchhoff stress
	tens3ds    m_T;				// LTE - 2nd Piola-Kirchhoff stress moment

	mat3ds     m_inf_str;		// LTE - infinitesimal strain
	tens3ds    m_inf_str_grad;	// LTE - infinitesimal strain gradient
	
	mat3ds     m_E;				// LTE - Green-Lagrange strain
	tens3ds    m_H;				// LTE - Green-Lagrange strain gradient
	
	mat3ds     m_e;				// LTE - Euler-Almansi strain
	tens3ds    m_h;				// LTE - Euler-Almansi strain graident

	double	   m_energy_diff;	// LTE - Difference between macro energy and volume averaged energy of RVE (should be zero) 

	tens4ds	   m_Ca;			//!< averaged material stiffness
	tens5ds    m_Da;
	tens6ds    m_Ea;
};

//-----------------------------------------------------------------------------
//! The micro-material implements material homogenization. The stress and tangents
//! are calculated by solving a micro-structural RVE problem and return the
//! averaged stress and condensed tangents.
//!
class FEMicroMaterial2O :	public FEElasticMaterial
{
public:
	FEMicroMaterial2O(FEModel* pfem);
	~FEMicroMaterial2O(void);

public:
	char	m_szrve[256];	//!< filename for RVE file
	char	m_szbc[256];	//!< name of nodeset defining boundary
	bool	m_bperiodic;	//!< periodic bc flag

protected:
	FEModel	m_rve;			//!< the master RVE (Representive Volume Element)
	bool	m_brve;			//!< flag indicating whether RVE was read in
	double	m_V0;			//!< initial volume of RVE
	vector<int> m_BN;		//!< boundary node flags

	double m_bb_x; double m_bb_y; double m_bb_z;  // LTE - RVE bounding box
	int m_num_ext_node;

public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
	void Stress2O(FEMaterialPoint &mp, int plot_on);

	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);
	void Tangent2O(FEMaterialPoint &mp, tens4ds& c, tens5ds& d, tens6ds& e);
	
	//! data initialization
	void Init();

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData();

protected:
	bool PrepRVE();
	bool PrepDisplacementBC();
	bool PrepPeriodicBC();
	void FindBoundaryNodes();

	void UpdateBC(FEModel& rve, mat3d& F, tens3drs& G);
	
	mat3ds AveragedStress(FEModel& rve, FEMaterialPoint &mp);
	void AveragedStress2O(FEModel& rve, FEMaterialPoint &mp, mat3ds &sa, tens3ds &taua);
	void AveragedStiffness(FEModel& rve, FEMaterialPoint &mp, tens4ds& c, tens5ds& d, tens6ds& e);

	void calculate_d2O(tens5ds& d, double K[3][3], double Ri[3], double Rj[3]);
	void calculate_e2O(tens6ds& e, double K[3][3], double Ri[3], double Rj[3]);
	
	void calc_energy_diff(FEModel& rve, FEMaterialPoint& pt);

	void AveragedStress2OPK1(FEModel& rve, FEMaterialPoint &mp, mat3d &PK1a, tens3drs &QK1a);
	void AveragedStress2OPK2(FEModel& rve, FEMaterialPoint &mp, mat3ds &Sa, tens3ds &Ta);

public:
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
