#pragma once
#include "FECore/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Multigenerational material point.
//! First generation exists at t=0. Second, third, etc. generations appear at t>0.
//! This material point stores the inverse of the relative deformation gradient of
//! second, third, etc. generations.  These relate the reference configuration of 
//! each generation relative to the first generation.

class FEMultigenerationMaterialPoint : public FEMaterialPoint
	{
	public:
		FEMultigenerationMaterialPoint(FEMaterialPoint* pt) : FEMaterialPoint(pt) {}
		
		FEMaterialPoint* Copy()
		{
			FEMultigenerationMaterialPoint* pt = new FEMultigenerationMaterialPoint(*this);
			if (m_pt) pt->m_pt = m_pt->Copy();
			return pt;
		}
		
		void Serialize(DumpFile& ar)
		{
			if (m_pt) m_pt->Serialize(ar);
			
			if (ar.IsSaving())
			{
				for (int i=0; i < (int) Fi.size(); i++)
					ar << Fi[i];
				for (int i=0; i < (int) Ji.size(); i++)
					ar << Ji[i];
			}
			else
			{
				for (int i=0; i < (int) Fi.size(); i++)
					ar >> Fi[i];
				for (int i=0; i < (int) Ji.size(); i++)
					ar >> Ji[i];
			}
			
			if (m_pt) m_pt->Serialize(ar);
		}
		
		void Init(bool bflag)
		{
			if (bflag)
			{
				Fi.clear();
				Ji.clear();
			}
			
			if (m_pt) m_pt->Init(bflag);
		}
		
	public:
		// multigenerational material data
		vector <mat3d> Fi;	//!< inverse of relative deformation gradient
		vector <double> Ji;	//!< determinant of Fi (store for efficiency)
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
		return new FEMultigenerationMaterialPoint(m_pMat[0]->CreateMaterialPointData());
	}
	
public:
	vector <FEElasticMaterial*>	m_pMat;	//!< pointers to elastic material
		
public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
		
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);
		
	//! data initialization and checking
	void Init();
		
	//! check existence of generation
	bool HasGeneration(const int igen);
	
	// declare the parameter list
//	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Global generation data

struct FEGenerationData {
	double	btime;	//!< generation birth time
	bool	born;	//!< flag for generation birth
};
