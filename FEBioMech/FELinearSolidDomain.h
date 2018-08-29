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
	FEMaterial* GetMaterial() override;

	//! Initialization
	bool Init() override;

	//! reset element data
	void Reset() override;

	//! initialize elements
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;

public: // overrides from FELinearElasticDomain

	//! Build the stiffness matrix
	void StiffnessMatrix(FESolver* psolver) override;

	// Calculate the RHS vector
	void RHS(FEGlobalVector& R) override;

	//! Update domain data
	void Update(const FETimeInfo& tp) override;

protected:
	void InitialStress(FESolidElement& el, vector<double>& fe);
	void InternalForce(FESolidElement& el, vector<double>& fe);

	void ElementStiffness(FESolidElement& el, matrix& ke);

protected:
	FESolidMaterial*	m_pMat;
};
