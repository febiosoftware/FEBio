#include "stdafx.h"
#include "FEPlotData.h"
#include "FEModel.h"
#include "tens3d.h"
#include "FESPRProjection.h"

//-----------------------------------------------------------------------------
FEPlotData::FEPlotData() : FECoreBase(FEPLOTDATA_ID)
{
	m_ntype = PLT_FLOAT;
	m_sfmt = FMT_NODE;
	m_nregion = FE_REGION_NODE;
	m_pfem = 0;

	m_arraySize = 0;
}

//-----------------------------------------------------------------------------
FEPlotData::FEPlotData(Region_Type R, Var_Type t, Storage_Fmt s) : FECoreBase(FEPLOTDATA_ID)
{ 
	m_ntype = t; 
	m_sfmt = s; 
    m_nregion = R;
	m_pfem = 0; 

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

//-----------------------------------------------------------------------------
void FEPlotData::Save(FEModel& fem, Archive& ar)
{
	switch (m_nregion)
	{
	case FE_REGION_NODE: SaveNodeData(fem, ar); break;
	case FE_REGION_DOMAIN: SaveDomainData(fem, ar); break;
	case FE_REGION_SURFACE: SaveSurfaceData(fem, ar); break;
	default:
		assert(false);
		break;
	}
}

//-----------------------------------------------------------------------------
void FEPlotData::SaveNodeData(FEModel &fem, Archive& ar)
{
	// store pointer to model
	m_pfem = &fem;

	// loop over all node sets
	// write now there is only one, namely the master node set
	// so we just pass the mesh
	int ndata = VarSize(DataType());

	int N = fem.GetMesh().Nodes();
	FEDataStream a; a.reserve(ndata*N);
	if (Save(fem.GetMesh(), a))
	{
		assert(a.size() == N*ndata);
		ar.WriteData(0, a.data());
	}
}

//-----------------------------------------------------------------------------
void FEPlotData::SaveDomainData(FEModel &fem, Archive& ar)
{
	// store pointer to model
	m_pfem = &fem;

	FEMesh& m = fem.GetMesh();
	int ND = m.Domains();

	// if the item list is empty, store all domains
	if (m_item.empty())
	{
		for (int i=0; i<ND; ++i) m_item.push_back(i);
	}

	// loop over all domains in the item list
	int N = (int)m_item.size();
	for (int i=0; i<ND; ++i)
	{
		// get the domain
		FEDomain& D = m.Domain(m_item[i]);

		// calculate the size of the data vector
		int nsize = VarSize(DataType());
		switch (m_sfmt)
		{
		case FMT_NODE: nsize *= D.Nodes(); break;
		case FMT_ITEM: nsize *= D.Elements(); break;
		case FMT_MULT:
			{
				// since all elements have the same type within a domain
				// we just grab the number of nodes of the first element 
				// to figure out how much storage we need
				FEElement& e = D.ElementRef(0);
				int n = e.Nodes();
				nsize *= n*D.Elements();
			}
			break;
		case FMT_REGION:
			// one value for this domain so nsize remains unchanged
			break;
		default:
			assert(false);
		}
		assert(nsize > 0);

		// fill data vector and save
		FEDataStream a; 
		a.reserve(nsize);
		if (Save(D, a))
		{
			assert(a.size() == nsize);
			ar.WriteData(m_item[i]+1, a.data());
		}
	}
}

//-----------------------------------------------------------------------------
//! Save surface data
//! \todo For the FMT_MULT option we are assuming 8 values per facet. I need to
//! make sure that the FEBioPlot assumes as many values.
void FEPlotData::SaveSurfaceData(FEModel &fem, Archive& ar)
{
	// store pointer to model
	m_pfem = &fem;

	// loop over all surfaces
	FEMesh& m = fem.GetMesh();
	int NS = m.Surfaces();
	for (int i=0; i<NS; ++i)
	{
		FESurface& S = m.Surface(i);

		// Determine data size.
		// Note that for the FMT_MULT case we are 
		// assuming 9 data entries per facet
		// regardless of the nr of nodes a facet really has
		// this is because for surfaces, all elements are not
		// necessarily of the same type
		// TODO: Fix the assumption of the FMT_MULT
		int nsize = VarSize(DataType());
		switch (m_sfmt)
		{
		case FMT_NODE: nsize *= S.Nodes(); break;
		case FMT_ITEM: nsize *= S.Elements(); break;
		case FMT_MULT: nsize *= 10*S.Elements(); break;
		case FMT_REGION: 
			// one value per surface so nsize remains unchanged
			break;
		default:
			assert(false);
		}

		// save data
		FEDataStream a; a.reserve(nsize);
		if (Save(S, a))
		{
			assert(a.size() == nsize);
			ar.WriteData(i+1, a.data());
		}
	}
}

//=================================================================================================

template <class T> void _writeNodalValuesT(FEDomain& dom, std::function<T(int)> f, FEDataStream& ar)
{
	int N = dom.Nodes();
	for (int i = 0; i<N; ++i) ar << f(i);
}

void writeNodalValues(FEDomain& dom, std::function<mat3ds(int)> f, FEDataStream& ar) { _writeNodalValuesT<mat3ds>(dom, f, ar); }

//=================================================================================================

template <class T> void _writeAverageElementValueT(FEDomain& dom, FEValuator<T>& var, FEDataStream& ar)
{
	// write solid element data
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
			s += var.eval(*el.GetMaterialPoint(j));
		}
		s *= f;

		ar << s;
	}
}

//-------------------------------------------------------------------------------------------------
void writeAverageElementValue(FEDomain& dom, FEValuator<double >& var, FEDataStream& ar) { _writeAverageElementValueT<double >(dom, var, ar); }
void writeAverageElementValue(FEDomain& dom, FEValuator<vec3d  >& var, FEDataStream& ar) { _writeAverageElementValueT<vec3d  >(dom, var, ar); }
void writeAverageElementValue(FEDomain& dom, FEValuator<mat3ds >& var, FEDataStream& ar) { _writeAverageElementValueT<mat3ds >(dom, var, ar); }
void writeAverageElementValue(FEDomain& dom, FEValuator<tens4ds>& var, FEDataStream& ar) { _writeAverageElementValueT<tens4ds>(dom, var, ar); }

//=================================================================================================

template <class Tin, class Tout> void _writeAverageElementValueT(FEDomain& dom, FEValuator<Tin>& var, FEDataStream& ar, std::function<Tout (const Tin& m)> flt)
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
			s += var.eval(*el.GetMaterialPoint(j));
		}
		s *= f;

		ar << flt(s);
	}
}

void writeAverageElementValue(FEDomain& dom, FEValuator<vec3d>& var, FEDataStream& ar, std::function<double(const vec3d& m)> flt)
{
	_writeAverageElementValueT<vec3d, double>(dom, var, ar, flt);
}

void writeAverageElementValue(FEDomain& dom, FEValuator<vec3d>& var, FEDataStream& ar, std::function<vec3d (const vec3d& m)> flt)
{
	_writeAverageElementValueT<vec3d, vec3d>(dom, var, ar, flt);
}

void writeAverageElementValue(FEDomain& dom, FEValuator<mat3ds>& var, FEDataStream& ar, std::function<double(const mat3ds& m)> flt)
{
	_writeAverageElementValueT<mat3ds, double>(dom, var, ar, flt);
}

void writeAverageElementValue(FEDomain& dom, FEValuator<mat3d>& var, FEDataStream& ar, std::function<double(const mat3d& m)> flt)
{
	_writeAverageElementValueT<mat3d, double>(dom, var, ar, flt);
}

void writeAverageElementValue(FEDomain& dom, FEValuator<tens3drs>& var, FEDataStream& ar, std::function<double(const tens3drs& m)> flt)
{
	_writeAverageElementValueT<tens3drs, double>(dom, var, ar, flt);
}

//=================================================================================================

template <class T> void _writeIntegratedElementValueT(FESolidDomain& dom, FEValuator<T>& var, FEDataStream& ar)
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
			T vj = var.eval(mp);
			double detJ = dom.detJ0(el, j)*gw[j];
			ew += vj*detJ;
		}
		ar << ew;
	}
}

//-------------------------------------------------------------------------------------------------
void writeIntegratedElementValue(FESolidDomain& dom, FEValuator<double>& var, FEDataStream& ar) { _writeIntegratedElementValueT<double>(dom, var, ar); }
void writeIntegratedElementValue(FESolidDomain& dom, FEValuator<vec3d>&  var, FEDataStream& ar) { _writeIntegratedElementValueT<vec3d >(dom, var, ar); }

//=================================================================================================
//-------------------------------------------------------------------------------------------------
void writeSPRElementValue(FESolidDomain& dom, FEValuator<mat3dd>& var, FEDataStream& ar, int interpolOrder)
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
			mat3dd v = var.eval(mp);

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
void writeSPRElementValue(FESolidDomain& dom, FEValuator<mat3ds>& var, FEDataStream& ar, int interpolOrder)
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
			mat3ds s = var.eval(mp);

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

//-------------------------------------------------------------------------------------------------
void writeNodalProjectedElementValues(FEDomain& dom, FEValuator<mat3ds>& var, FEDataStream& ar)
{
	// tensor component look-up table
	int LUT[6][2] = { { 0,0 },{ 1,1 },{ 2,2 },{ 0,1 },{ 1,2 },{ 0,2 } };

	// temp storage 
	mat3ds s[FEElement::MAX_NODES];
	double si[FEElement::MAX_INTPOINTS];
	double sn[FEElement::MAX_NODES];

	// loop over all elements
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int ne = e.Nodes();
		int ni = e.GaussPoints();

		// loop over tensor components
		for (int j = 0; j<6; ++j)
		{
			// get the integration point values
			int j0 = LUT[j][0];
			int j1 = LUT[j][1];
			for (int k = 0; k<ni; ++k)
			{
				FEMaterialPoint& mp = *e.GetMaterialPoint(k);
				mat3ds s = var.eval(mp);
				si[k] = s(j0, j1);
			}

			// project to nodes
			e.project_to_nodes(si, sn);

			// store stress component
			for (int k = 0; k<ne; ++k) s[k](j0, j1) = sn[k];
		}

		// push data to archive
		for (int j = 0; j<ne; ++j) ar << s[j];
	}
}
