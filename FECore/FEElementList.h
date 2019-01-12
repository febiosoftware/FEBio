#pragma once
#include "fecore_api.h"

class FEMesh;
class FEElement;

//-----------------------------------------------------------------------------
//! utitlity class for accessing all elements without having to go throug the domains
class FEElementList
{
public:
	class iterator
	{
	public:
		iterator() { m_pmesh = 0; m_ndom = -1; m_nel = -1; }
		iterator(FEMesh* pm) { m_pmesh = pm; m_ndom = 0; m_nel = 0; }

		FECORE_API FEElement& operator*();

		FECORE_API void operator ++ ();

		bool operator != (iterator& it)
		{
			return ((m_ndom!=it.m_ndom)||(m_nel != it.m_nel));
		}

	public:
		FEMesh*		m_pmesh;	// pointer to mesh
		int			m_ndom;		// domain index
		int			m_nel;		// element index
	};

public:
	FEElementList(FEMesh& m) : m_mesh(m){}

	iterator begin() { return iterator(&m_mesh); }
	iterator end() { return iterator(); }

private:
	FEMesh&	m_mesh;
};
