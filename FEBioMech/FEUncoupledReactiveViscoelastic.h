#pragma once
#include "FEUncoupledMaterial.h"
#include "FEBondRelaxation.h"
#include "FEReactiveVEMaterialPoint.h"

//-----------------------------------------------------------------------------
//! This class implements a large deformation reactive viscoelastic material
//! with uncoupled strain energy density formulation
//
class FEUncoupledReactiveViscoelasticMaterial :	public FEUncoupledMaterial
{
public:
    //! default constructor
    FEUncoupledReactiveViscoelasticMaterial(FEModel* pfem);
        
    //! get the elastic base material
	FEUncoupledMaterial* GetBaseMaterial() { return m_pBase; }
    
    //! Set the base material
    void SetBaseMaterial(FEUncoupledMaterial* pbase) { m_pBase = pbase; }
    
    //! get the elastic bond material
	FEUncoupledMaterial* GetBondMaterial() { return m_pBond; }
    
    //! Set the base material
    void SetBondMaterial(FEUncoupledMaterial* pbond) { m_pBond = pbond; }
    
public:
    //! data initialization
    bool Init() override;
    
    //! stress function
    mat3ds DevStress(FEMaterialPoint& pt) override;
    
    //! tangent function
    tens4ds DevTangent(FEMaterialPoint& pt) override;
    
    //! strain energy density function
    double DevStrainEnergyDensity(FEMaterialPoint& pt) override;
    
    //! cull generations
    void CullGenerations(FEMaterialPoint& pt);
    
    //! evaluate bond mass fraction for a given generation
    double BreakingBondMassFraction(FEMaterialPoint& pt, const int ig, const mat3ds D);
    
    //! evaluate bond mass fraction of reforming generation
    double ReformingBondMassFraction(FEMaterialPoint& pt);
    
    //! detect new generation
    bool NewGeneration(FEMaterialPoint& pt);
    
    //! returns a pointer to a new material point object
    FEMaterialPoint* CreateMaterialPointData() override;
    
private:
    FEUncoupledMaterial*	m_pBase;	//!< pointer to elastic solid material for polymer base
	FEUncoupledMaterial*	m_pBond;	//!< pointer to elastic solid material for reactive bonds
	FEBondRelaxation*		m_pRelx;    //!< pointer to bond relaxation material for reactive bonds
    
public:
    double	m_wmin;		//!< minimum value of relaxation
    int     m_btype;    //!< bond kinetics type
    int     m_ttype;    //!< bond breaking trigger type
    
    DECLARE_FECORE_CLASS();
};
