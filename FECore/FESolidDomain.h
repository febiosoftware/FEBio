#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! abstract base class for 3D volumetric elements
class FESolidDomain : public FEDomain
{
public:
	//! constructor
	FESolidDomain(int ntype, FEMesh* pm, FEMaterial* pmat) : FEDomain(ntype, pm, pmat) {}

	//! create storage for elements
	void create(int nsize) { m_Elem.resize(nsize); }

	//! return nr of elements
	int Elements() { return m_Elem.size(); }

	//! element access
	FESolidElement& Element(int n) { return m_Elem[n]; }
	FESolidElement& ElementRef(int n) { return m_Elem[n]; }

	int GetElementType() { return m_Elem[0].Type(); }

	bool Initialize(FEModel& fem);

	int Nodes() { return (int) m_Node.size(); }
	FENode& Node(int i);

	//! find the element in which point y lies
	FESolidElement* FindElement(vec3d y, double r[3]);

	//! Calculate deformation gradient at integration point n
	double defgrad(FESolidElement& el, mat3d& F, int n);

protected:
	vector<int>				m_Node;		//!< node list
	vector<FESolidElement>	m_Elem;		//!< array of elements
};
