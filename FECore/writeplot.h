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
#include "FEMesh.h"
#include "FESurface.h"
#include "FEDataStream.h"
#include "FESolidDomain.h"
#include "FEDomainParameter.h"
#include "fecore_api.h"
#include <functional>

//=================================================================================================
template <class T> void writeNodalValues(FEMesh& mesh, FEDataStream& ar, std::function<T(const FENode& node)> f)
{
	for (int i = 0; i<mesh.Nodes(); ++i) ar << f(mesh.Node(i));
}

//=================================================================================================
template <class T> void writeNodalValues(FEMeshPartition& dom, FEDataStream& ar, std::function<T(int)> f)
{
	for (int i = 0; i<dom.Nodes(); ++i) ar << f(i);
}

//=================================================================================================
template <class T> void writeElementValue(FEMeshPartition& dom, FEDataStream& ar, std::function<T(int nface)> f)
{
	for (int i = 0; i<dom.Elements(); ++i) ar << f(i);
}

//=================================================================================================
template <class T> void writeNodalValues(FENodeSet& nset, FEDataStream& ar, std::function<T (const FEMaterialPoint&)> var)
{
	FEMesh& mesh = *nset.GetMesh();
	vector<T> data(mesh.Nodes(), T(0.0));

	FEMaterialPoint mp;
	for (int i = 0; i < nset.Size(); ++i) {
		FENode& node = mesh.Node(nset[i]);
		mp.m_r0 = node.m_r0;
		mp.m_index = i;
		data[nset[i]] = var(mp);
	}
	ar << data;
}

//=================================================================================================
template <class T> void writeIntegratedElementValue(FESurface& surf, FEDataStream& ar, std::function<T(const FEMaterialPoint& mp)> fnc)
{
	T s(0.0);
	for (int i = 0; i<surf.Elements(); ++i) {
		FESurfaceElement& el = surf.Element(i);
		double* w = el.GaussWeights();
		for (int j = 0; j<el.GaussPoints(); ++j) s += fnc(*el.GetMaterialPoint(j))*w[j];
	}
	ar << s;
}

//=================================================================================================
template <class T> void writeElementValue(FEMeshPartition& dom, FEDataStream& ar, std::function<T(const FEMaterialPoint& mp)> fnc)
{
	int NE = dom.Elements();
	std::vector<T> v(NE);
#pragma omp parallel for shared(v)
	for (int i = 0; i<NE; ++i) {
		FEElement& el = dom.ElementRef(i);
		v[i] = fnc(*el.GetMaterialPoint(0));
	}

	for (int i = 0; i < NE; ++i)
		ar << v[i];
}

//=================================================================================================
template <class T> void writeAverageElementValue(FEMeshPartition& dom, FEDataStream& ar, std::function<T(const FEMaterialPoint& mp)> fnc)
{
	int NE = dom.Elements();
	std::vector<T> v(NE);
#pragma omp parallel for shared(v)
	for (int i = 0; i<NE; ++i) {
		FEElement& el = dom.ElementRef(i);
		T s(0.0);
		for (int j = 0; j<el.GaussPoints(); ++j) s += fnc(*el.GetMaterialPoint(j));
		v[i] = s / (double)el.GaussPoints();
	}

	for (int i=0; i<NE; ++i)
		ar << v[i];
}

//=================================================================================================
template <class T> void writeAverageElementValue(FEMeshPartition& dom, FEDataStream& ar, std::function<T(FEElement& el, int ip)> fnc)
{
	for (int i = 0; i<dom.Elements(); ++i) {
		FEElement& el = dom.ElementRef(i);
		T s(0.0);
		for (int j = 0; j<el.GaussPoints(); ++j) s += fnc(el, j);
		ar << s / (double) el.GaussPoints();
	}
}

//=================================================================================================
template <class Tin, class Tout> void writeAverageElementValue(FEMeshPartition& dom, FEDataStream& ar, std::function<Tin(const FEMaterialPoint&)> fnc, std::function<Tout(const Tin& m)> flt)
{
	for (int i = 0; i<dom.Elements(); ++i) {
		FEElement& el = dom.ElementRef(i);
		Tin s(0.0);
		for (int j = 0; j<el.GaussPoints(); ++j) s += fnc(*el.GetMaterialPoint(j));
		ar << flt(s / (double) el.GaussPoints());
	}
}

//=================================================================================================
template <class Tin, class Tout> void writeAverageElementValue(FEMeshPartition& dom, FEDataStream& ar, std::function<Tin(FEElement& el, int ip)> fnc, std::function<Tout(const Tin& m)> flt)
{
	for (int i = 0; i<dom.Elements(); ++i) {
		FEElement& el = dom.ElementRef(i);
		Tin s(0.0);
		for (int j = 0; j<el.GaussPoints(); ++j) s += fnc(el, j);
		ar << flt(s / (double)el.GaussPoints());
	}
}

//=================================================================================================
template <class T> void writeAverageElementValue(FEMeshPartition& dom, FEDataStream& ar, FEDomainParameter* var)
{
	for (int i = 0; i<dom.Elements(); ++i) {
		FEElement& el = dom.ElementRef(i);
		T s(0.0);
		for (int j = 0; j < el.GaussPoints(); ++j)
		{
			FEParamValue v = var->value(*el.GetMaterialPoint(j));
			s += v.value<T>();
		}
		ar << s / (double)el.GaussPoints();
	}
}

//=================================================================================================
template <class T> void writeIntegratedElementValue(FESolidDomain& dom, FEDataStream& ar, std::function<T(const FEMaterialPoint& mp)> fnc)
{
	for (int i = 0; i<dom.Elements(); ++i) {
		FESolidElement& el = dom.Element(i);
		double* gw = el.GaussWeights();

		T ew(0.0);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			ew += fnc(mp)*dom.detJ0(el, j)*gw[j];
		}
		ar << ew;
	}
}

//=================================================================================================
template <class T> void writeNodalProjectedElementValues(FEMeshPartition& dom, FEDataStream& ar, std::function<T(const FEMaterialPoint&)> var)
{
	// temp storage 
	T si[FEElement::MAX_INTPOINTS];
	T sn[FEElement::MAX_NODES];

	// loop over all elements
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int ne = e.Nodes();
		int ni = e.GaussPoints();

		// get the integration point values
		for (int k = 0; k<ni; ++k)
		{
			FEMaterialPoint& mp = *e.GetMaterialPoint(k);
			T s = var(mp);
			si[k] = s;
		}

		// project to nodes
		e.project_to_nodes(si, sn);

		// push data to archive
		for (int j = 0; j<ne; ++j) ar << sn[j];
	}
}

//=================================================================================================
template <class T> void writeNodalProjectedElementValues(FESurface& dom, FEDataStream& ar, std::function<T(const FEMaterialPoint&)> var)
{
	T gi[FEElement::MAX_INTPOINTS];
	T gn[FEElement::MAX_NODES];

	// loop over all the elements in the domain
	int NE = dom.Elements();
	for (int i = 0; i < NE; ++i)
	{
		// get the element and loop over its integration points
		// we only calculate the element's average
		// but since most material parameters can only defined 
		// at the element level, this should get the same answer
		FESurfaceElement& e = dom.Element(i);
		int nint = e.GaussPoints();
		int neln = e.Nodes();

		for (int j = 0; j < nint; ++j)
		{
			// get the material point data for this integration point
			FEMaterialPoint& mp = *e.GetMaterialPoint(j);
			gi[j] = var(mp);
		}

		e.FEElement::project_to_nodes(gi, gn);

		// store the result
		for (int j = 0; j < neln; ++j) ar << gn[j];
	}
}

//-----------------------------------------------------------------------------
// helper functions for writing SPR projected element values
// TODO: I needed to give these functions a different name because of the implicit conversion between mat3ds and mat3dd
FECORE_API void writeSPRElementValueMat3dd(FESolidDomain& dom, FEDataStream& ar, std::function<mat3dd(const FEMaterialPoint&)> fnc, int interpolOrder = -1);
FECORE_API void writeSPRElementValueMat3ds(FESolidDomain& dom, FEDataStream& ar, std::function<mat3ds(const FEMaterialPoint&)> fnc, int interpolOrder = -1);

// Helper functions for mapping data
FECORE_API void ProjectToNodes(FEDomain& dom, vector<double>& nodeVals, function<double(FEMaterialPoint& mp)> f);
FECORE_API void writeRelativeError(FEDomain& dom, FEDataStream& a, function<double(FEMaterialPoint& mp)> f);
