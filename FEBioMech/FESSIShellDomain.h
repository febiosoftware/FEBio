#pragma once
#include <FECore/FEShellDomain.h>

//-----------------------------------------------------------------------------
// This class extends the FEShellDomain and implements the solid-shell interface (SSI) logic.
// It is used by the new shell formulation.
class FESSIShellDomain : public FEShellDomain
{
public:
	FESSIShellDomain(FEMesh* mesh);

	//! Update element data prior to solving time step
	void PreSolveUpdate(const FETimeInfo& timeInfo);

protected:
	//! Find interfaces between solid element faces and shell elements
	void FindSSI();

public:
	//! calculates covariant basis vectors at an integration point
	void CoBaseVectors0(FEShellElement& el, int n, vec3d g[3]);

	//! calculates contravariant basis vectors at an integration point
	void ContraBaseVectors0(FEShellElement& el, int n, vec3d g[3]);

	// inverse jacobian with respect to reference frame
	double invjac0(FEShellElement& el, double J[3][3], int n);

	// jacobian with respect to reference frame
	double detJ0(FEShellElement& el, int n);

protected:
	bool	m_binit;
};
