#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! Abstract base class for shell elements
class FEShellDomain : public FEDomain
{
public:
	//! constructor
	FEShellDomain(FEMesh* pm) : FEDomain(FE_DOMAIN_SHELL, pm) {}

	//! create storage for elements
	void Create(int nsize, int elemType);

	//! return nr of elements
	int Elements() const { return (int)m_Elem.size(); }

	//! element access
	FEShellElement& Element(int n) { return m_Elem[n]; }
	FEElement& ElementRef(int n) { return m_Elem[n]; }

	int GetElementType() { return m_Elem[0].Type(); }

	//! Update element data prior to solving time step
	void PreSolveUpdate(const FETimeInfo& timeInfo);

	//! Reset element data
	void Reset();

    //! calculates covariant basis vectors at an integration point
    void CoBaseVectors0(FEShellElement& el, int n, vec3d g[3]);
    
    //! calculates contravariant basis vectors at an integration point
    void ContraBaseVectors0(FEShellElement& el, int n, vec3d g[3]);
    
	// inverse jacobian with respect to reference frame
	double invjac0(FEShellElement& el, double J[3][3], int n);

	// jacobian with respect to reference frame
	double detJ0(FEShellElement& el, int n);

	//! Serialize domain data to archive
	void Serialize(DumpStream& ar);

public:
    //! Find interfaces between solid element faces and shell elements
    void FindSSI();
    
protected:
	vector<FEShellElement>	m_Elem;	//!< array of elements
    bool                    m_binit;    //!< initialization flag
};
