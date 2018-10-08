#include "stdafx.h"
#include "FEPlotData.h"
#include "FEModel.h"
#include "tens3d.h"
#include "FESPRProjection.h"
#include "FESurface.h"
#include "FEDomain.h"
#include "FESolidDomain.h"

//-----------------------------------------------------------------------------
FEPlotData::FEPlotData(FEModel* fem) : FECoreBase(FEPLOTDATA_ID)
{
	m_ntype = PLT_FLOAT;
	m_sfmt = FMT_NODE;
	m_nregion = FE_REGION_NODE;
	m_fem = fem;

	m_arraySize = 0;
}

//-----------------------------------------------------------------------------
FEPlotData::FEPlotData(FEModel* fem, Region_Type R, Var_Type t, Storage_Fmt s) : FECoreBase(FEPLOTDATA_ID)
{ 
	m_ntype = t; 
	m_sfmt = s; 
    m_nregion = R;
	m_fem = fem; 

	m_arraySize = 0;
}

//-----------------------------------------------------------------------------
void FEPlotData::SetArraySize(int n)
{
	m_arraySize = n;
}

//-----------------------------------------------------------------------------
int FEPlotData::GetArraysize() const
{
	return m_arraySize;
}

//-----------------------------------------------------------------------------
int FEPlotData::VarSize(Var_Type t)
{
	int ndata = 0;
	switch (DataType())
	{
	case PLT_FLOAT  : ndata =  1; break;
	case PLT_VEC3F  : ndata =  3; break;
	case PLT_MAT3FS : ndata =  6; break;
	case PLT_MAT3FD : ndata =  3; break;
    case PLT_TENS4FS: ndata = 21; break;
	case PLT_MAT3F  : ndata =  9; break;
	case PLT_ARRAY  : ndata = GetArraysize(); break;
	case PLT_ARRAY_VEC3F: ndata = GetArraysize()*3; break;
	}
	assert(ndata);
	return ndata;
}

//-----------------------------------------------------------------------------
void FEPlotData::SetDomainName(const char* szdom)
{
	strcpy(m_szdom, szdom); 
}

//=================================================================================================

template <class T> void _writeNodalValues(FEMesh& mesh, FEDataStream& ar, std::function<T (const FENode& node)> f)
{
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		ar << f(node);
	}
}

void writeNodalValues(FEMesh& mesh, FEDataStream& ar, std::function<double (const FENode& node)> f) { _writeNodalValues<double>(mesh, ar, f); }
void writeNodalValues(FEMesh& mesh, FEDataStream& ar, std::function<vec3d  (const FENode& node)> f) { _writeNodalValues<vec3d >(mesh, ar, f); }

//=================================================================================================

template <class T> void _writeNodalValuesT(FEDomain& dom, FEDataStream& ar, std::function<T(int)> f)
{
	int N = dom.Nodes();
	for (int i = 0; i<N; ++i) ar << f(i);
}

void writeNodalValues(FEDomain& dom, FEDataStream& ar, std::function<double(int)> fnc) { _writeNodalValuesT<double>(dom, ar, fnc); }
void writeNodalValues(FEDomain& dom, FEDataStream& ar, std::function<mat3ds(int)> fnc) { _writeNodalValuesT<mat3ds>(dom, ar, fnc); }

//=================================================================================================
template <class T> void _writeElementValue(FEDomain& dom, FEDataStream& ar, std::function<T (int nface)> f)
{
	int NF = dom.Elements();
	for (int i = 0; i<NF; ++i)
	{
		ar << f(i);
	}
}

void writeElementValue(FEDomain& dom, FEDataStream& ar, std::function<double(int nface)> f) { _writeElementValue<double>(dom, ar, f); }
void writeElementValue(FEDomain& dom, FEDataStream& ar, std::function<vec3d (int nface)> f) { _writeElementValue<vec3d >(dom, ar, f); }

//=================================================================================================
template <class T> void _writeAverageElementValueT(FEDomain& dom, FEDataStream& ar, std::function<T (const FEMaterialPoint& mp)> fnc)
{
	int N = dom.Elements();
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		T s(0.0);
		int nint = el.GaussPoints();
		double f = 1.0 / (double)nint;

		// we output the average value values of the gauss points
		for (int j = 0; j<nint; ++j)
		{
			s += fnc(*el.GetMaterialPoint(j));
		}
		s *= f;

		ar << s;
	}
}

//-------------------------------------------------------------------------------------------------
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<double  (const FEMaterialPoint& mp)> fnc) { _writeAverageElementValueT<double >(dom, ar, fnc); }
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<vec3d   (const FEMaterialPoint& mp)> fnc) { _writeAverageElementValueT<vec3d  >(dom, ar, fnc); }
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<mat3ds  (const FEMaterialPoint& mp)> fnc) { _writeAverageElementValueT<mat3ds >(dom, ar, fnc); }
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<tens4ds (const FEMaterialPoint& mp)> fnc) { _writeAverageElementValueT<tens4ds>(dom, ar, fnc); }


//=================================================================================================
template <class T> void _writeAverageElementValueT(FEDomain& dom, FEDataStream& ar, std::function<T (FEElement& el, int ip)> fnc)
{
	int N = dom.Elements();
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		T s(0.0);
		int nint = el.GaussPoints();
		double f = 1.0 / (double)nint;

		// we output the average value values of the gauss points
		for (int j = 0; j<nint; ++j)
		{
			s += fnc(el, j);
		}
		s *= f;

		ar << s;
	}
}

//-----------------------------------------------------------------------------
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<mat3ds(FEElement& el, int ip)> fnc) { _writeAverageElementValueT(dom, ar, fnc); }

//=================================================================================================
template <class T> void _writeSummedElementValue(FEDomain& dom, FEDataStream& ar, std::function<T (const FEMaterialPoint& mp)> fnc)
{
	int N = dom.Elements();
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		T s(0.0);
		int nint = el.GaussPoints();

		// we output the sum of the integration point values
		for (int j = 0; j<nint; ++j)
		{
			s += fnc(*el.GetMaterialPoint(j));
		}

		ar << s;
	}
}

//-------------------------------------------------------------------------------------------------
void writeSummedElementValue(FEDomain& dom, FEDataStream& ar, std::function<double(const FEMaterialPoint& mp)> fnc)
{
	_writeSummedElementValue<double>(dom, ar, fnc);
}

void writeSummedElementValue(FEDomain& dom, FEDataStream& ar, std::function<vec3d (const FEMaterialPoint& mp)> fnc)
{
	_writeSummedElementValue<vec3d>(dom, ar, fnc);
}


//=================================================================================================

template <class Tin, class Tout> void _writeAverageElementValueT(FEDomain& dom, FEDataStream& ar, std::function<Tin (const FEMaterialPoint&)> fnc, std::function<Tout (const Tin& m)> flt)
{
	// write solid element data
	int N = dom.Elements();
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		Tin s(0.0);
		int nint = el.GaussPoints();
		double f = 1.0 / (double)nint;

		// we output the average value values of the gauss points
		for (int j = 0; j<nint; ++j)
		{
			s += fnc(*el.GetMaterialPoint(j));
		}
		s *= f;

		ar << flt(s);
	}
}

void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<vec3d (const FEMaterialPoint&)> fnc, std::function<double(const vec3d& m)> flt)
{
	_writeAverageElementValueT<vec3d, double>(dom, ar, fnc, flt);
}

void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<vec3d(const FEMaterialPoint&)> fnc, std::function<vec3d (const vec3d& m)> flt)
{
	_writeAverageElementValueT<vec3d, vec3d>(dom, ar, fnc, flt);
}

void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<mat3d (const FEMaterialPoint&)> fnc, std::function<double(const mat3d& m)> flt)
{
	_writeAverageElementValueT<mat3d, double>(dom, ar, fnc, flt);
}

void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<mat3ds (const FEMaterialPoint&)> fnc, std::function<double(const mat3ds& m)> flt)
{
	_writeAverageElementValueT<mat3ds, double>(dom, ar, fnc, flt);
}

void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<tens3drs (const FEMaterialPoint&)> fnc, std::function<double(const tens3drs& m)> flt)
{
	_writeAverageElementValueT<tens3drs, double>(dom, ar, fnc, flt);
}


//=================================================================================================
template <class Tin, class Tout> void _writeAverageElementValueT(FEDomain& dom, FEDataStream& ar, std::function<Tin(FEElement& el, int ip)> fnc, std::function<Tout(const Tin& m)> flt)
{
	// write solid element data
	int N = dom.Elements();
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		Tin s(0.0);
		int nint = el.GaussPoints();
		double f = 1.0 / (double)nint;

		// we output the average value values of the gauss points
		for (int j = 0; j<nint; ++j)
		{
			s += fnc(el, j);
		}
		s *= f;

		ar << flt(s);
	}
}


//-------------------------------------------------------------------------------------------------
void writeAverageElementValue(FEDomain& dom, FEDataStream& ar, std::function<mat3ds(FEElement& el, int ip)> fnc, std::function<double(const mat3ds& m)> flt)
{
	_writeAverageElementValueT<mat3ds, double>(dom, ar, fnc, flt);
}

//=================================================================================================

template <class T> void _writeIntegratedElementValueT(FESolidDomain& dom, FEDataStream& ar, std::function<T (const FEMaterialPoint& mp)> fnc)
{
	for (int i = 0; i<dom.Elements(); ++i)
	{
		FESolidElement& el = dom.Element(i);
		double* gw = el.GaussWeights();

		// integrate
		T ew(0.0);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			T vj = fnc(mp);
			double detJ = dom.detJ0(el, j)*gw[j];
			ew += vj*detJ;
		}
		ar << ew;
	}
}

//-------------------------------------------------------------------------------------------------
void writeIntegratedElementValue(FESolidDomain& dom, FEDataStream& ar , std::function<double (const FEMaterialPoint& mp)> fnc) { _writeIntegratedElementValueT<double>(dom, ar, fnc); }
void writeIntegratedElementValue(FESolidDomain& dom, FEDataStream& ar , std::function<vec3d  (const FEMaterialPoint& mp)> fnc) { _writeIntegratedElementValueT<vec3d >(dom, ar, fnc); }

//=================================================================================================
//-------------------------------------------------------------------------------------------------
void writeSPRElementValueMat3dd(FESolidDomain& dom, FEDataStream& ar, std::function<mat3dd(const FEMaterialPoint&)> fnc, int interpolOrder)
{
	int NN = dom.Nodes();
	int NE = dom.Elements();

	// build the element data array
	vector< vector<double> > ED[3];
	ED[0].resize(NE);
	ED[1].resize(NE);
	ED[2].resize(NE);
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& e = dom.Element(i);
		int nint = e.GaussPoints();
		ED[0][i].assign(nint, 0.0);
		ED[1][i].assign(nint, 0.0);
		ED[2][i].assign(nint, 0.0);
	}

	// this array will store the results
	FESPRProjection map;
	map.SetInterpolationOrder(interpolOrder);
	vector<double> val[3];

	// fill the ED array
	for (int i = 0; i < NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		for (int j = 0; j < nint; ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			mat3dd v = fnc(mp);

			ED[0][i][j] = v.diag(0);
			ED[1][i][j] = v.diag(1);
			ED[2][i][j] = v.diag(2);
		}
	}

	// project to nodes
	map.Project(dom, ED[0], val[0]);
	map.Project(dom, ED[1], val[1]);
	map.Project(dom, ED[2], val[2]);

	// copy results to archive
	for (int i = 0; i<NN; ++i)
	{
		ar.push_back((float)val[0][i]);
		ar.push_back((float)val[1][i]);
		ar.push_back((float)val[2][i]);
	}
}

//-------------------------------------------------------------------------------------------------
void writeSPRElementValueMat3ds(FESolidDomain& dom, FEDataStream& ar, std::function<mat3ds(const FEMaterialPoint&)> fnc, int interpolOrder)
{
	const int LUT[6][2] = { { 0,0 },{ 1,1 },{ 2,2 },{ 0,1 },{ 1,2 },{ 0,2 } };

	int NN = dom.Nodes();
	int NE = dom.Elements();

	// build the element data array
	vector< vector<double> > ED[6];
	for (int n = 0; n < 6; ++n)
	{
		ED[n].resize(NE);
		for (int i = 0; i < NE; ++i)
		{
			FESolidElement& e = dom.Element(i);
			int nint = e.GaussPoints();
			ED[n][i].assign(nint, 0.0);
		}
	}

	// this array will store the results
	FESPRProjection map;
	map.SetInterpolationOrder(interpolOrder);
	vector<double> val[6];

	// fill the ED array
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			mat3ds s = fnc(mp);

			// loop over stress components
			for (int n = 0; n < 6; ++n)
			{
				ED[n][i][j] = s(LUT[n][0], LUT[n][1]);
			}
		}
	}

	// project to nodes
	// loop over stress components
	for (int n = 0; n<6; ++n)
	{
		map.Project(dom, ED[n], val[n]);
	}

	// copy results to archive
	for (int i = 0; i<NN; ++i)
	{
		ar.push_back((float)val[0][i]);
		ar.push_back((float)val[1][i]);
		ar.push_back((float)val[2][i]);
		ar.push_back((float)val[3][i]);
		ar.push_back((float)val[4][i]);
		ar.push_back((float)val[5][i]);
	}
}

//=================================================================================================
template <class T> void _writeNodalProjectedElementValues(FEDomain& dom, FEDataStream& ar, std::function<T (const FEMaterialPoint&)> var)
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

//-------------------------------------------------------------------------------------------------
void writeNodalProjectedElementValues(FEDomain& dom, FEDataStream& ar, std::function<double(const FEMaterialPoint&)> fnc)
{
	_writeNodalProjectedElementValues<double>(dom, ar, fnc);
}

//-------------------------------------------------------------------------------------------------
void writeNodalProjectedElementValues(FEDomain& dom, FEDataStream& ar, std::function<vec3d (const FEMaterialPoint&)> fnc)
{
	_writeNodalProjectedElementValues<vec3d>(dom, ar, fnc);
}

//-------------------------------------------------------------------------------------------------
void writeNodalProjectedElementValues(FEDomain& dom, FEDataStream& ar, std::function<mat3ds(const FEMaterialPoint&)> fnc)
{
	_writeNodalProjectedElementValues<mat3ds>(dom, ar, fnc);
}

//=================================================================================================
template <class T> void _writeNodalProjectedElementValues(FESurface& dom, FEDataStream& ar, std::function<T (const FEMaterialPoint&)> var)
{
	T gi[FEElement::MAX_INTPOINTS];
	T gn[FEElement::MAX_NODES    ];

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
		// NOTE: Note that we always need to store 10 entries. This is because of a limitation of the plot file format.
		for (int j = 0; j < 10; ++j) ar << gn[j];
	}
}

//-------------------------------------------------------------------------------------------------
void writeNodalProjectedElementValues(FESurface& dom, FEDataStream& ar, std::function<double(const FEMaterialPoint&)> var) 
{
	_writeNodalProjectedElementValues<double>(dom, ar, var);
}

void writeNodalProjectedElementValues(FESurface& dom, FEDataStream& ar, std::function<vec3d (const FEMaterialPoint&)> var)
{
	_writeNodalProjectedElementValues<vec3d>(dom, ar, var);
}

//=================================================================================================
void writeNodalValues(FENodeSet& nset, FEDataStream& ar, std::function<double(const FEMaterialPoint&)> var)
{
	FEMesh& mesh = *nset.GetMesh();
	int NN = mesh.Nodes();
	vector<double> data(NN, 0.0);

	const std::vector<int> nodeList = nset.GetNodeList();
	FEMaterialPoint mp;
	for (int i = 0; i < nset.size(); ++i)
	{
		int nodeId = nodeList[i];
		FENode& node = mesh.Node(nodeId);
		mp.m_r0 = node.m_r0;
		mp.m_index = i;

		double vi = var(mp);

		data[nodeId] = vi;
	}
	ar << data;
}
