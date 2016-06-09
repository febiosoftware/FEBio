#pragma once
#include "FECore/FESolidDomain.h"

//-----------------------------------------------------------------------------
class FESolidMaterial;

//-----------------------------------------------------------------------------
class FELinearElasticDomain
{
public:
	FELinearElasticDomain(FEModel* pfem);
	virtual ~FELinearElasticDomain(){}
	virtual void StiffnessMatrix(FESolver* psolver) = 0;
	virtual void RHS(FEGlobalVector& R) = 0;
};

//-----------------------------------------------------------------------------
//! Class describing a linear elastic solid domain
class FELinearSolidDomain : public FESolidDomain, public FELinearElasticDomain
{
public:
	//! constructor
	FELinearSolidDomain(FEModel* pfem, FEMaterial* pmat);

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial();

	//! Initialization
	bool Initialize(FEModel& fem);

	//! reset element data
	void Reset();

	//! initialize elements
	void PreSolveUpdate(const FETimeInfo& timeInfo);

public: // overrides from FELinearElasticDomain

	//! Build the stiffness matrix
	void StiffnessMatrix(FESolver* psolver);

	// Calculate the RHS vector
	void RHS(FEGlobalVector& R);

	//! Update domain data
	void Update(const FETimeInfo& tp);

protected:
	void InitialStress(FESolidElement& el, vector<double>& fe);
	void InternalForce(FESolidElement& el, vector<double>& fe);

	void ElementStiffness(FESolidElement& el, matrix& ke);

protected:
	FESolidMaterial*	m_pMat;
};
