/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include <FECore/FEShellDomain.h>
#include <FECore/FEModelParam.h>
#include <functional>
#include "febiomech_api.h"
#include <FECore/FEDofList.h>
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

	//! serialization
	void Serialize(DumpStream& ar) override;
    
	//! Update element data prior to solving time step
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
    //! Initialize shell normals
    void InitShells() override;
    
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
    
    //! calculates covariant basis vectors at any point
    void CoBaseVectors(FEShellElement& el, double r, double s, double t, vec3d g[3], const double alpha);
    
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
    
    //! calculates contravariant basis vectors at any point
    void ContraBaseVectors(FEShellElement& el, double r, double s, double t, vec3d g[3], const double alpha);
    
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
	FEDofList	m_dofU;		// displacement dofs
	FEDofList	m_dofSU;	// shell displacement dofs
	FEDofList	m_dofR;		// rigid rotation
    
public:
    bool        m_bnodalnormals; // flag for nodal (true) or element (false) normals

	DECLARE_FECORE_CLASS();
};


//-----------------------------------------------------------------------------
void writeIntegratedElementValue(FESSIShellDomain& dom, FEDataStream& ar, std::function<double (const FEMaterialPoint& mp)> fnc);
void writeIntegratedElementValue(FESSIShellDomain& dom, FEDataStream& ar, std::function<vec3d  (const FEMaterialPoint& mp)> fnc);
