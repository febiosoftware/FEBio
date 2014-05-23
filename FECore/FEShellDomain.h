#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! Abstract base class for shell elements
class FEShellDomain : public FEDomain
{
public:
	//! constructor
	FEShellDomain(int ntype, FEMesh* pm) : FEDomain(ntype, FE_DOMAIN_SHELL, pm) {}

	//! create storage for elements
	void create(int nsize) { m_Elem.resize(nsize); }

	//! return nr of elements
	int Elements() { return m_Elem.size(); }

	//! element access
	FEShellElement& Element(int n) { return m_Elem[n]; }
	FEElement& ElementRef(int n) { return m_Elem[n]; }

	int GetElementType() { return m_Elem[0].Type(); }

	//! Initialize elements
	void InitElements();

	bool Initialize(FEModel& fem);

	//! Reset element data
	void Reset();

	int Nodes() { return (int) m_Node.size(); }
	FENode& Node(int i);

	// calculate deformation gradient
	double defgrad(FEShellElement& el, mat3d& F, int n);

	// inverse jacobian with respect to reference frame
	double invjac0(FEShellElement& el, double J[3][3], int n);

	// inverse jacobian with respect to current frame
	double invjact(FEShellElement& el, double J[3][3], int n);

	// jacobian with respect to reference frame
	double detJ0(FEShellElement& el, int n);

public:
	//! shallow copy
	void ShallowCopy(DumpStream& dmp, bool bsave);

	//! Serialize domain data to archive
	void Serialize(DumpFile& ar);

protected:
	vector<int>				m_Node;	//!< node list
	vector<FEShellElement>	m_Elem;	//!< array of elements
};
