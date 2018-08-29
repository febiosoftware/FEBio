#pragma once
#include "FEMultiphasic.h"

//-----------------------------------------------------------------------------
// forward declaration of material class
class FEMultiphasicMultigeneration;

//-----------------------------------------------------------------------------
//! Multigenerational SBM material point.
//! This material point stores the inverse of the relative deformation gradient,
//! and the increment in the mass of solid-bound molecular species in
//! multiple generations.

class FEMultigenSBMMaterialPoint : public FEMaterialPoint
{
public:
	FEMultigenSBMMaterialPoint(FEMultiphasicMultigeneration* pm, FEMaterialPoint* pt) : m_pmat(pm), FEMaterialPoint(pt) { m_tgen = 0.0; }
    
	FEMaterialPoint* Copy();
    
	void Serialize(DumpStream& ar);
    
	void Init();
	void Update(const FETimeInfo& timeInfo);
    
public:
	// multigenerational material data
    int                         m_ngen;     //!< number of generations
	vector <mat3d>              m_Fi;       //!< inverse of relative deformation gradient
	vector <double>             m_Ji;       //!< determinant of Fi (store for efficiency)
    int                         m_nsbm;     //!< number of solid-bound molecules
    vector< vector<double> >    m_gsbmr;    //!< sbmr content at each generation
    vector< vector<double> >    m_gsbmrp;   //!< gsbmr at previous time point
    vector<double>              m_lsbmr;    //!< last generation sbmr values
	double                      m_tgen;     //!< last generation time
    
private:
	FEMultiphasicMultigeneration*	m_pmat;
};

//-----------------------------------------------------------------------------
//! Multigeneration multiphasic material.

class FEMultiphasicMultigeneration : public FEMultiphasic
{
public:
	//! constructor
	FEMultiphasicMultigeneration(FEModel* pfem);
    
    //! returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() override;
	
    //! Update solid bound molecules
    void UpdateSolidBoundMolecules(FEMaterialPoint& mp) override;

	int CheckGeneration(const double t);
    double GetGenerationTime(const int igen);
    
public:
	double	m_gtime;	//!< time duration of each generation
    
    DECLARE_PARAMETER_LIST();

};
