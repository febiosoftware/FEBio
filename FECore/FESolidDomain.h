#pragma once
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! abstract base class for 3D volumetric elements
class FESolidDomain : public FEDomain
{
public:
	//! constructor
	FESolidDomain(FEMesh* pm) : FEDomain(FE_DOMAIN_SOLID, pm) {}

	//! create storage for elements
	void create(int nsize) { m_Elem.resize(nsize); }

	//! return nr of elements
	int Elements() { return m_Elem.size(); }

	//! reset data (overridden from FEDomain)
	void Reset();

	//! element access
	FESolidElement& Element(int n) { return m_Elem[n]; }
	FEElement& ElementRef(int n) { return m_Elem[n]; }

	int GetElementType() { return m_Elem[0].Type(); }

	bool Initialize(FEModel& fem);

	int Nodes() { return (int) m_Node.size(); }
	FENode& Node(int i);
	int NodeIndex(int i) { return m_Node[i]; }

	//! find the element in which point y lies
	FESolidElement* FindElement(vec3d y, double r[3]);

	//! Calculate deformation gradient at integration point n
	double defgrad(FESolidElement& el, mat3d& F, int n);

	//! Calculate deformation gradient at integration point n
	double defgrad(FESolidElement& el, mat3d& F, double r, double s, double t);

	//! calculate inverse jacobian matrix w.r.t. reference frame
	double invjac0(FESolidElement& el, double J[3][3], int n);

	//! calculate inverse jacobian matrix w.r.t. reference frame
	double invjac0(FESolidElement& el, double J[3][3], double r, double s, double t);

	//! calculate inverse jacobian matrix w.r.t. current frame
	double invjact(FESolidElement& el, double J[3][3], int n);

	//! calculate inverse jacobian matrix w.r.t. reference frame
	double invjact(FESolidElement& el, double J[3][3], double r, double s, double t);

	//! calculate gradient of function at integration points
	vec3d gradient(FESolidElement& el, double* fn, int n);

	//! calculate gradient of function at integration points
	vec3d gradient(FESolidElement& el, vector<double>& fn, int n);

	//! calculate jacobian in reference frame
	double detJ0(FESolidElement& el, int n);

	//! calculate jacobian in current frame
	double detJt(FESolidElement& el, int n);

public:
	//! serialize data to archive
	void Serialize(DumpFile& ar);

	//! shallow copy
	void ShallowCopy(DumpStream& dmp, bool bsave);

protected:
	vector<int>				m_Node;		//!< node list
	vector<FESolidElement>	m_Elem;		//!< array of elements
};
