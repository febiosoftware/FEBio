// FENodeElemList.cpp: implementation of the FENodeElemList class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FENodeElemList.h"
#include "FESurface.h"
#include "FEMesh.h"
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! This function builds the node-element list for a surface

void FENodeElemList::Create(FESurface& s)
{
	int i, j, n;

	// get the number of nodes
	int nn = s.Nodes();

	// get the number of elements
	int ne = s.Elements();

	// create nodal valence array
	m_nval.resize(nn);
	m_pn.resize(nn);

	// clear valence table
	m_nval.zero();

	// fill valence table
	int nsize = 0;
	for (i=0; i<ne; ++i)
	{
		FESurfaceElement& el = s.Element(i);

		for (j=0; j<el.Nodes(); ++j)
		{
			n = el.m_lnode[j];
			m_nval[n]++;
			nsize++;
		}
	}

	// create the element reference array
	m_eref.resize(nsize);

	// set eref pointers
	m_pn[0] = 0;
	for (i=1; i<nn; ++i)
	{
		m_pn[i] = m_pn[i-1] + m_nval[i-1];
	}

	// reset valence pointers
	for (i=0; i<nn; ++i) m_nval[i] = 0;

	// fill eref table
	for (i=0; i<ne; ++i)
	{
		FESurfaceElement& el = s.Element(i);

		for (j=0; j<el.Nodes(); ++j)
		{
			n = el.m_lnode[j];
			m_eref[m_pn[n] + m_nval[n]] = &el;
			m_nval[n]++;
		}
	}
}

//-----------------------------------------------------------------------------
//! This function builds the node-element list for a mesh

void FENodeElemList::Create(FEMesh& mesh)
{
	int i, j, n, nd;

	// get the number of nodes
	int NN = mesh.Nodes();

	// create nodal valence array
	m_nval.assign(NN, 0);
	m_pn.resize(NN);

	// fill valence table
	int nsize = 0;
	for (nd=0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& d = mesh.Domain(nd);
		for (i=0; i<d.Elements(); ++i)
		{
			FEElement& el = d.ElementRef(i);
			for (j=0; j<el.Nodes(); ++j)
			{
				n = el.m_node[j];
				m_nval[n]++;
				nsize++;
			}
		}
	}

	// create the element reference array
	m_eref.resize(nsize);

	// set eref pointers
	m_pn[0] = 0;
	for (i=1; i<NN; ++i)
	{
		m_pn[i] = m_pn[i-1] + m_nval[i-1];
	}

	// reset valence pointers
	for (i=0; i<NN; ++i) m_nval[i] = 0;

	// fill eref table
	for (nd=0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& d = mesh.Domain(nd);
		for (i=0; i<d.Elements(); ++i)
		{
			FEElement& el = d.ElementRef(i);
			for (j=0; j<el.Nodes(); ++j)
			{
				n = el.m_node[j];
				m_eref[m_pn[n] + m_nval[n]] = &el;
				m_nval[n]++;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function builds the node-element list for a domain

void FENodeElemList::Create(FEDomain& dom)
{
	int i, j, n;

	// get the mesh
	FEMesh& mesh = *dom.GetMesh();

	// get the number of nodes
	int NN = mesh.Nodes();

	// create nodal valence array
	m_nval.assign(NN, 0);
	m_pn.resize(NN);

	// fill valence table
	int nsize = 0;
	for (i=0; i<dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		for (j=0; j<el.Nodes(); ++j)
		{
			n = el.m_node[j];
			m_nval[n]++;
			nsize++;
		}
	}

	// create the element reference array
	m_eref.resize(nsize);

	// set eref pointers
	m_pn[0] = 0;
	for (i=1; i<NN; ++i)
	{
		m_pn[i] = m_pn[i-1] + m_nval[i-1];
	}

	// reset valence pointers
	for (i=0; i<NN; ++i) m_nval[i] = 0;

	// fill eref table
	for (i=0; i<dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);
		for (j=0; j<el.Nodes(); ++j)
		{
			n = el.m_node[j];
			m_eref[m_pn[n] + m_nval[n]] = &el;
			m_nval[n]++;
		}
	}
}
