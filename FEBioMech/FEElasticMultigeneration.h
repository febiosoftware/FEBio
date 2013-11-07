#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Global generation data

struct FEGenerationData {
	double	btime;	//!< generation birth time
	bool	born;	//!< flag for generation birth
};

//-----------------------------------------------------------------------------
// forward declaration of material class
class FEElasticMultigeneration;

//-----------------------------------------------------------------------------
//! Multigenerational material point.
//! First generation exists at t=0. Second, third, etc. generations appear at t>0.
//! This material point stores the inverse of the relative deformation gradient of
//! second, third, etc. generations.  These relate the reference configuration of 
//! each generation relative to the first generation.

class FEMultigenerationMaterialPoint : public FEMaterialPoint
{
public:
	FEMultigenerationMaterialPoint(FEElasticMultigeneration* pm, FEMaterialPoint* pt) : m_pmat(pm), FEMaterialPoint(pt) { m_tgen = 0.0; }
		
	FEMaterialPoint* Copy();
		
	void Serialize(DumpFile& ar);
		
	void Init(bool bflag);
		
public:
	// multigenerational material data
	vector <mat3d> Fi;	//!< inverse of relative deformation gradient
	vector <double> Ji;	//!< determinant of Fi (store for efficiency)
	double	m_tgen;		//!< last generation time
private:
	FEElasticMultigeneration*	m_pmat;
};

//-----------------------------------------------------------------------------
//! Multigenerational solid

class FEElasticMultigeneration : public FEElasticMaterial
{
public:
	FEElasticMultigeneration () {}
		
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() 
	{ 
		// use the zero-th generation material point as the base elastic material point
		return new FEMultigenerationMaterialPoint(this, m_pMat[0]->CreateMaterialPointData());
	}

	void AddMaterial(FEElasticMaterial* pmat);
	
public:
	vector <FEElasticMaterial*>	m_pMat;	//!< pointers to elastic material

public:
	//! return material properties
	int Properties() { return (int) m_pMat.size(); }

	//! return a material property
	FEMaterial* GetProperty(int i) { return m_pMat[i]; }
		
public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
		
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);
		
	//! data initialization and checking
	void Init();
		
	//! check existence of generation
	bool HasGeneration(const int igen);

	static void PushGeneration(FEGenerationData* G);
	static int CheckGeneration(const double t);

public:
	static vector<FEGenerationData*> m_MG;	//!< multigeneration data

	// declare the parameter list
//	DECLARE_PARAMETER_LIST();
};
