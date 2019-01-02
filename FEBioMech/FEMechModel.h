#pragma once
#include <FECore/FEModel.h>

//---------------------------------------------------------------------------------------
class FERigidSystem;

//---------------------------------------------------------------------------------------
// This class extends the basic FEModel class by adding a rigid body system
class FEMechModel : public FEModel
{
public:
	FEMechModel();

	// get the rigid system
	FERigidSystem* GetRigidSystem();

	// clear all model data
	void Clear() override;

	// initialize the rigid system
	bool InitRigidSystem() override;

	// model activation
	void Activate() override;

	// reset
	bool Reset() override;

	//! Initialize shells
	void InitShells() override;

	// find a parameter value
	FEParamValue GetParameterValue(const ParamString& param) override;

	//! find a model componnet from its class ID
	FEModelComponent* FindModelComponent(int nid) override;

	//! serialize data for restarts
	void Serialize(DumpStream& ar) override;

	//! Build the matrix profile for this model
	void BuildMatrixProfile(FEGlobalMatrix& G, bool breset) override;

private:
	FERigidSystem*	m_prs;

	DECLARE_FECORE_CLASS();
};
