#pragma once
#include <FECore/FEShellDomain.h>
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
// This class extends the FEShellDomain and implements the solid-shell interface (SSI) logic.
// It is used by the new shell formulation.
class FESSIShellDomain : public FEShellDomain
{
public:
	FESSIShellDomain(FEModel* pfem);

    //! initialize domain
    //! one-time initialization, called during model initialization
    bool Initialize();
    
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

public:
    //! calculates covariant basis vectors at an integration point
    void CoBaseVectors(FEShellElement& el, int n, vec3d g[3]);
    
    //! calculates contravariant basis vectors at an integration point
    void ContraBaseVectors(FEShellElement& el, int n, vec3d g[3]);
    
    // jacobian with respect to current frame
    double detJ(FEShellElement& el, int n);
    
    // calculate deformation gradient
    double defgrad(FEShellElement& el, mat3d& F, int n);
    
    // inverse jacobian with respect to current frame
    double invjact(FEShellElement& el, double J[3][3], int n);
    
    //! evaluate a vector function over the shell
    vec3d evaluate(FEShellElement& el, vec3d* vn, vec3d* dvn, int n);
    
    //! evaluate a scalar function over the shell
    double evaluate(FEShellElement& el, double* pn, double* dpn, int n);
    
    //! calculate the gradient of a scalar function over the shell
    vec3d gradient(FEShellElement& el, double* pn, double* dpn, int n);
    
protected:
    int     m_dofx;
    int     m_dofy;
    int     m_dofz;
    int     m_dofu;
    int     m_dofv;
    int     m_dofw;
};
