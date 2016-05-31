#pragma once
#include "FECore/FEMaterial.h"
#include "FEBioMech/FEBodyForce.h"
#include "FEElasticFluid.h"
#include "FEViscousFluid.h"

//-----------------------------------------------------------------------------
//! Fluid material point class.
//
class FEFluidMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
    FEFluidMaterialPoint();

	//! create a shallow copy
	FEMaterialPoint* Copy();

	//! data serialization
	void Serialize(DumpStream& ar);

	//! Data initialization
	void Init(bool bflag);

public:
    mat3ds RateOfDeformation() { return m_L.sym(); }
    mat3da Spin() { return m_L.skew(); }
    vec3d  Vorticity() { return vec3d(m_L(2,1)-m_L(1,2), m_L(0,2)-m_L(2,0), m_L(1,0)-m_L(0,1)); }
    
public:
    // fluid data
    vec3d       m_r0;       //!< material position
    vec3d       m_vt;       //!< velocity
    vec3d       m_vp;       //!< velocity at previous time
    vec3d       m_at;       //!< acceleration
    mat3d       m_L;        //!< velocity gradient
    double      m_J;        //!< determinant of fluid deformation gradient
    double      m_Jp;       //!< determinant of fluid deformation gradient at previous time
    vec3d       m_gradJ;    //!< gradient of J
	double		m_p;		//!< elastic fluid pressure
    mat3ds		m_s;		//!< fluid stress
    vec3d       m_lapv;     //!< laplacian of velocity
    vec3d       m_gdiv;     //!< gradient of divergence of velocity
};

//-----------------------------------------------------------------------------
//! Base class for fluid materials.

class FEFluid : public FEMaterial
{
public:
	FEFluid(FEModel* pfem);
	
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();
	
public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
	
    //! tangent of stress with respect to strain J
    mat3ds Tangent_Strain(FEMaterialPoint& mp);
    
    //! tangent of stress with respect to rate of deformation tensor D
    tens4ds Tangent_RateOfDeformation(FEMaterialPoint& mp)  { return m_pViscous->Tangent_RateOfDeformation(mp); }
    
	//! referential fluid density
	double ReferentialDensity() { return m_rhor; }

    //! calculate current fluid density
    double Density(FEMaterialPoint& pt);
    
    //! return elastic part
    FEElasticFluid* GetElastic() { return m_pElastic; }
    
    //! return viscous part
    FEViscousFluid* GetViscous() { return m_pViscous; }
    
    //! kinematic viscosity
    double KinematicViscosity(FEMaterialPoint& mp);
    
    //! acoustic speed
    double AcousticSpeed(FEMaterialPoint& mp);
    
private: // material properties
	FEPropertyT<FEElasticFluid>	m_pElastic;	//!< pointer to elastic part fluid material
    FEPropertyT<FEViscousFluid> m_pViscous; //!< pointer to viscous part of fluid material
	
public:
    double						m_rhor;     //!< referential fluid density
    bool                        m_bsupg;    //!< flag for SUPG stabilization
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};
