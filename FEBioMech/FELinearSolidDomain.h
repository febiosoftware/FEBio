#pragma once
#include "FECore/FESolidDomain.h"
#include "FECore/FESolver.h"
#include "FECore/FEModel.h"
#include "FESolidMaterial.h"

//-----------------------------------------------------------------------------
class FELinearElasticDomain
{
public:
	virtual ~FELinearElasticDomain(){}
	virtual void StiffnessMatrix(FESolver* psolver) = 0;
	virtual void RHS(FEGlobalVector& R) = 0;
	virtual void UpdateStresses(FEModel& fem) = 0;
};

//-----------------------------------------------------------------------------
//! Class describing a linear elastic solid domain
class FELinearSolidDomain : public FESolidDomain, public FELinearElasticDomain
{
public:
	//! constructor
	FELinearSolidDomain(FEMesh* pm, FEMaterial* pmat);

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() { return m_pMat; }

	//! Initialization
	bool Initialize(FEModel& fem);

	//! activate
	void Activate();

	//! Unpack solid element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! reset element data
	void Reset();

	//! initialize elements
	void InitElements();

public: // overrides from FELinearElasticDomain

	//! Build the stiffness matrix
	void StiffnessMatrix(FESolver* psolver);

	// Calculate the RHS vector
	void RHS(FEGlobalVector& R);

	//! Update the element stresses
	void UpdateStresses(FEModel& fem);

protected:
	void InitialStress(FESolidElement& el, vector<double>& fe);
	void InternalForce(FESolidElement& el, vector<double>& fe);

	void ElementStiffness(FESolidElement& el, matrix& ke);

protected:
	FESolidMaterial*	m_pMat;
};
