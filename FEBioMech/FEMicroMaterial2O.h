#pragma once
#include "FEElasticMaterial2O.h"
#include <FECore/tens3d.h>
#include <FECore/tens4d.h>
#include <FECore/tens5d.h>
#include <FECore/tens6d.h>
#include "FE2OMicroConstraint.h"
#include "FEMicroMaterial.h"
#include "FERVEModel2O.h"

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
	mat3d	m_Pa;				//!< averaged PK1 stress

	tens3drs   m_G;				//!< gradient of deformation gradient
	tens3drs   m_Qa;			//!< average higher-order stress tensor

	tens4d	m_Ca;				//!< averaged 4-th order stiffness
	tens5d	m_La, m_Ha;			//!< averaged 5-th order stiffness
	tens6d	m_Ja;				//!< averaged 6-th order stiffness

/*	tens3ds    m_tau;			// TODO: remove

	mat3ds     m_S;				// LTE - 2nd Piola-Kirchhoff stress (TODO: remove)
	tens3ds    m_T;				// LTE - 2nd Piola-Kirchhoff stress moment (TODO: remove)

	tens3ds    m_inf_str_grad;	// LTE - infinitesimal strain gradient (TODO: remove)
	
	mat3ds     m_E;				// LTE - Green-Lagrange strain	(TODO: remove)
	tens3ds    m_H;				// LTE - Green-Lagrange strain gradient (TODO: remove)
	
	mat3ds     m_e;				// LTE - Euler-Almansi strain (TODO: remove)
	tens3ds    m_h;				// LTE - Euler-Almansi strain graident (TODO: remove)

	double     m_macro_energy;	// LTE - Total macroscopic strain energy (TODO: remove)
	double	   m_micro_energy;	// LTE - Total volume-average of strain energy throughout the RVE solution (TODO: remove)
	double	   m_energy_diff;	// LTE - Difference between macro energy and volume averaged energy of RVE (should be zero)  (TODO: remove)

	double	   m_macro_energy_inc;	// LTE - Macroscopic strain energy increment (TODO: remove)
	double	   m_micro_energy_inc;	// LTE - Microscopic strain energy increment (TODO: remove)

	tens4ds	   m_Ca;			//!< Averaged rank 4 material stiffness
	tens5ds    m_Da;			//!< Averaged rank 5 material stiffness
	tens6ds    m_Ea;			//!< Averaged rank 6 material stiffness

	mat3d		m_F_prev;		//!< deformation gradient at previous converged time step
	tens3drs	m_G_prev; 		//!< gradient of deformation gradient at previous time point
*/

	FEMicroModel2O m_rve;				//!< local copy of the rve		
};

//-----------------------------------------------------------------------------
//! The micro-material implements material homogenization. The stress and tangents
//! are calculated by solving a micro-structural RVE problem and return the
//! averaged stress and condensed tangents.
//!
class FEMicroMaterial2O : public FEElasticMaterial2O
{
public:
	FEMicroMaterial2O(FEModel* pfem);
	~FEMicroMaterial2O(void);

public:
	char			m_szrve[256];	//!< filename for RVE file
	char			m_szbc[256];	//!< name of nodeset defining boundary
	bool			m_bperiodic;	//!< periodic bc flag
	FERVEModel2O	m_mrve;			//!< the master RVE (Representive Volume Element)

public:
	//! calculate stress at material point
	void Stress2O(FEMaterialPoint &mp);

	//! calculate tangent stiffness at material point
	void Tangent2O(FEMaterialPoint &mp, tens4d& C, tens5d& L, tens5d& H, tens6d& J);
	
	//! data initialization
	bool Init();

	//! create material point data
	FEMaterialPoint* CreateMaterialPointData();

protected:
	void calc_energy_diff(FEModel& rve, FEMaterialPoint& pt);

public:
	int Probes() { return (int) m_probe.size(); }
	FEMicroProbe& Probe(int i) { return *m_probe[i]; }

private:
	// We don't need these functions
	// TODO: Perhaps we should derive this class directly from FEMaterial so
	//       that these functions are not inherited. 
	mat3ds Stress(FEMaterialPoint& pt);
	tens4ds Tangent(FEMaterialPoint& pt);

protected:
	FEVecPropertyT<FEMicroProbe>	m_probe;

public:
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
