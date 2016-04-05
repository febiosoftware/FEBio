#pragma once
#include "FEElasticMaterial2O.h"
#include "FECore/FEModel.h"
#include "FECore/FEMaterial.h"
#include "FEBioMech/FEPeriodicBoundary2O.h"
#include "FECore/tens3d.h"
#include "FECore/tens4d.h"
#include "FECore/tens5d.h"
#include "FECore/tens6d.h"
#include "FE2OMicroConstraint.h"
#include "FEMicroMaterial.h"
#include "FERVEModel.h"

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
	void Serialize(DumpStream& ar);

public:
	tens3drs   m_G;				// LTE - Deformation Hessian
	tens3ds    m_tau;			// LTE - Cauchy stress moment

	tens3drs   m_G_prev; 		// LTE - Deformation Hessian

	mat3d      m_PK1;			// LTE - 1st Piola-Kirchhoff stress
	tens3drs   m_QK1;			// LTE - 1st Piola-Kirchhoff stress moment

	mat3ds     m_S;				// LTE - 2nd Piola-Kirchhoff stress
	tens3ds    m_T;				// LTE - 2nd Piola-Kirchhoff stress moment

	tens3ds    m_inf_str_grad;	// LTE - infinitesimal strain gradient
	
	mat3ds     m_E;				// LTE - Green-Lagrange strain
	tens3ds    m_H;				// LTE - Green-Lagrange strain gradient
	
	mat3ds     m_e;				// LTE - Euler-Almansi strain
	tens3ds    m_h;				// LTE - Euler-Almansi strain graident

	double     m_macro_energy;	// LTE - Total macroscopic strain energy
	double	   m_micro_energy;	// LTE - Total volume-average of strain energy throughout the RVE solution
	double	   m_energy_diff;	// LTE - Difference between macro energy and volume averaged energy of RVE (should be zero) 

	double	   m_macro_energy_inc;	// LTE - Macroscopic strain energy increment
	double	   m_micro_energy_inc;	// LTE - Microscopic strain energy increment

	bool m_rve_init;			// LTE - Flag indicating that the rve has been initialized
	FEModel m_rve;				// LTE - Current copy of the rve		
	FEModel m_rve_prev;			// LTE - Previous converged state of the rve

	tens4ds	   m_Ca;			// LTE - Averaged rank 4 material stiffness
	tens5ds    m_Da;			// LTE - Averaged rank 5 material stiffness
	tens6ds    m_Ea;			// LTE - Averaged rank 6 material stiffness
};

//-----------------------------------------------------------------------------
//! The micro-material implements material homogenization. The stress and tangents
//! are calculated by solving a micro-structural RVE problem and return the
//! averaged stress and condensed tangents.
//!
class FEMicroMaterial2O :	public FEElasticMaterial2O
{
public:
	FEMicroMaterial2O(FEModel* pfem);
	~FEMicroMaterial2O(void);

public:
	char		m_szrve[256];	//!< filename for RVE file
	char		m_szbc[256];	//!< name of nodeset defining boundary
	bool		m_bperiodic;	//!< periodic bc flag
	FERVEModel	m_mrve;			//!< the master RVE (Representive Volume Element)

public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
	void Stress2O(FEMaterialPoint &mp);

	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);
	void Tangent2O(FEMaterialPoint &mp, tens4ds& c, tens5ds& d, tens6ds& e);
	
	//! data initialization
	bool Init();

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData();

protected:
	void UpdateBC(FEModel& rve, mat3d& F, tens3drs& G);
	
	mat3ds AveragedStress(FEModel& rve, FEMaterialPoint &mp);
	void AveragedStress2O(FEModel& rve, FEMaterialPoint &mp, mat3ds &sa, tens3ds &taua);
	void AveragedStiffness(FEModel& rve, FEMaterialPoint &mp, tens4ds& c, tens5ds& d, tens6ds& e);

	void calculate_d2O(tens5ds& d, double K[3][3], double Ri[3], double Rj[3]);
	double calc_5ds_comp(double K[3][3], double Ri[3], double Rj[3], int i, int j, int k, int l, int m);

	void calculate_e2O(tens6ds& e, double K[3][3], double Ri[3], double Rj[3]);
	double calc_6ds_comp(double K[3][3], double Ri[3], double Rj[3], int i, int j, int k, int l, int m, int n);

	void calc_energy_diff(FEModel& rve, FEMaterialPoint& pt);

	void AveragedStress2OPK1(FEModel& rve, FEMaterialPoint &mp, mat3d &PK1a, tens3drs &QK1a);
	void AveragedStress2OPK2(FEModel& rve, FEMaterialPoint &mp, mat3ds &Sa, tens3ds &Ta);

public:
	int Probes() { return (int) m_probe.size(); }
	FEMicroProbe& Probe(int i) { return *m_probe[i]; }

protected:
	FEVecPropertyT<FEMicroProbe>	m_probe;

public:
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
