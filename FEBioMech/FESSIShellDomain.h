#pragma once
#include <FECore/FEShellDomain.h>
#include <FECore/FEModel.h>
#include <FECore/FEModelParam.h>
#include <functional>
#include "febiomech_api.h"
class FEDataStream;

//-----------------------------------------------------------------------------
// This class extends the FEShellDomain and implements the solid-shell interface (SSI) logic.
// It is used by the new shell formulation.
class FEBIOMECH_API FESSIShellDomain : public FEShellDomainNew
{
public:
	FESSIShellDomain(FEModel* pfem);

    //! initialize domain
    //! one-time initialization, called during model initialization
	bool Init() override;
    
	//! Update element data prior to solving time step
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
protected:
	//! Find interfaces between solid element faces and shell elements
	void FindSSI();

public:
	//! calculates covariant basis vectors at an integration point
	void CoBaseVectors0(FEShellElement& el, int n, vec3d g[3]);

    //! calculates covariant basis vectors at any point
    void CoBaseVectors0(FEShellElement& el, double r, double s, double t, vec3d g[3]);

	//! calculates contravariant basis vectors at an integration point
	void ContraBaseVectors0(FEShellElement& el, int n, vec3d g[3]);

    //! calculates contravariant basis vectors at any point
    void ContraBaseVectors0(FEShellElement& el, double r, double s, double t, vec3d g[3]);

	// inverse jacobian with respect to reference frame at an integration point
	double invjac0(FEShellElement& el, double J[3][3], int n);

    // inverse jacobian with respect to reference frame at any point
    double invjac0(FEShellElement& el, double J[3][3], double r, double s, double t);
    
	// jacobian with respect to reference frame
	double detJ0(FEShellElement& el, int n);

    // jacobian with respect to current frame at any point
    double detJ0(FEShellElement& el, double r, double s, double t);
    
public:
    //! calculates covariant basis vectors at an integration point
    void CoBaseVectors(FEShellElement& el, int n, vec3d g[3]);
    
    //! calculates covariant basis vectors at any point
    void CoBaseVectors(FEShellElement& el, double r, double s, double t, vec3d g[3]);
    
    //! calculates covariant basis vectors at an integration point at previous time
    void CoBaseVectorsP(FEShellElement& el, int n, vec3d g[3]);
    
    //! calculates covariant basis vectors at an integration point at intermediate time
    void CoBaseVectors(FEShellElement& el, int n, vec3d g[3], const double alpha);
    
    //! calculates contravariant basis vectors at an integration point
    void ContraBaseVectors(FEShellElement& el, int n, vec3d g[3]);
    
    //! calculates contravariant basis vectors at an integration point at intermediate time
    void ContraBaseVectors(FEShellElement& el, int n, vec3d g[3], const double alpha);
    
    //! calculates contravariant basis vectors at any point
    void ContraBaseVectors(FEShellElement& el, double r, double s, double t, vec3d g[3]);
    
    // jacobian with respect to current frame at an integration point
    double detJ(FEShellElement& el, int n);
    
    // jacobian with respect to current frame at an integration point at intermediate time
    double detJ(FEShellElement& el, int n, const double alpha);
    
    // jacobian with respect to current frame at any point
    double detJ(FEShellElement& el, double r, double s, double t);
    
    // calculate deformation gradient at an integration point
    double defgrad(FEShellElement& el, mat3d& F, int n);
    
    // calculate deformation gradient at any point
    double defgrad(FEShellElement& el, mat3d& F, double r, double s, double t);
    
    // calculate deformation gradient at an integration point at previous time
    double defgradp(FEShellElement& el, mat3d& F, int n);
    
    // inverse jacobian with respect to current frame
    double invjact(FEShellElement& el, double J[3][3], int n);
    
    //! evaluate a vector function over the shell
    vec3d evaluate(FEShellElement& el, vec3d* vn, vec3d* dvn, int n);
    
    //! evaluate a scalar function over the shell
    double evaluate(FEShellElement& el, double* pn, double* dpn, int n);
    
    //! calculate the gradient of a scalar function over the shell
    vec3d gradient(FEShellElement& el, double* pn, double* dpn, int n);
    
    //! evaluate a scalar function over the shell
    double evaluate(FEShellElement& el, vector<double> pn, vector<double> dpn, int n);
    
    //! calculate the gradient of a scalar function over the shell
    vec3d gradient(FEShellElement& el, vector<double> pn, vector<double> dpn, int n);
    
	//! Functions for element-DOF updates
	virtual void UpdateEAS(vector<double>& ui) {}
	virtual void UpdateIncrementsEAS(vector<double>& ui, const bool binc) {}

	void Update(const FETimeInfo& tp) override;

protected:
    int     m_dofx;
    int     m_dofy;
    int     m_dofz;
    int     m_dofsx;
    int     m_dofsy;
    int     m_dofsz;
    int     m_dofsxp;
    int     m_dofsyp;
    int     m_dofszp;
};


//-----------------------------------------------------------------------------
void writeIntegratedElementValue(FESSIShellDomain& dom, FEDataStream& ar, std::function<double (const FEMaterialPoint& mp)> fnc);
void writeIntegratedElementValue(FESSIShellDomain& dom, FEDataStream& ar, std::function<vec3d  (const FEMaterialPoint& mp)> fnc);
