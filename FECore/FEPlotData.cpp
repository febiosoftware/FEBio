#include "stdafx.h"
#include "FEPlotData.h"

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
