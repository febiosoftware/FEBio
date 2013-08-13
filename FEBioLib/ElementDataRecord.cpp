#include "stdafx.h"
#include "ElementDataRecord.h"
#include "FEBioMix/FETriphasic.h"
#include "FEBioMix/FEMultiphasic.h"

//-----------------------------------------------------------------------------
double FELogElemPosX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_rt.x;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemPosY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_rt.y;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemPosZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_rt.z;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemJacobian::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_J;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.xx();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.yy();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.zz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainXY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.xy();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainYZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.yz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrainXZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		val += E.xz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrain1::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		E.exact_eigen(l);
		val += l[0];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrain2::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		E.exact_eigen(l);
		val += l[1];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStrain3::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		mat3ds E = pt.Strain();
		E.exact_eigen(l);
		val += l[2];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.xx();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.yy();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.zz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressXY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.xy();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressYZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.yz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStressXZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_s.xz();
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStress1::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		pt.m_s.exact_eigen(l);
		val += l[0];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStress2::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		pt.m_s.exact_eigen(l);
		val += l[1];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemStress3::value(FEElement& el)
{
	double l[3];
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		pt.m_s.exact_eigen(l);
		val += l[2];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientXX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(0,0);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientXY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(0,1);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientXZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(0,2);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientYX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(1,0);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientYY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(1,1);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientYZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(1,2);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientZX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(2,0);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientZY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(2,1);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemDeformationGradientZZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEElasticMaterialPoint& pt = *el.m_State[i]->ExtractData<FEElasticMaterialPoint>();
		val += pt.m_F(2,2);
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemFluidPressure::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEBiphasicMaterialPoint* ppt = el.m_State[i]->ExtractData<FEBiphasicMaterialPoint>();
		if (ppt) val += ppt->m_pa;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemFluidFluxX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEBiphasicMaterialPoint* ppt = el.m_State[i]->ExtractData<FEBiphasicMaterialPoint>();
		if (ppt) val += ppt->m_w.x;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemFluidFluxY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEBiphasicMaterialPoint* ppt = el.m_State[i]->ExtractData<FEBiphasicMaterialPoint>();
		if (ppt) val += ppt->m_w.y;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemFluidFluxZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEBiphasicMaterialPoint* ppt = el.m_State[i]->ExtractData<FEBiphasicMaterialPoint>();
		if (ppt) val += ppt->m_w.z;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteConcentration::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESoluteMaterialPoint* ppt = el.m_State[i]->ExtractData<FESoluteMaterialPoint>();
		if (ppt) val += ppt->m_ca;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteFluxX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESoluteMaterialPoint* ppt = el.m_State[i]->ExtractData<FESoluteMaterialPoint>();
		if (ppt) val += ppt->m_j.x;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteFluxY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESoluteMaterialPoint* ppt = el.m_State[i]->ExtractData<FESoluteMaterialPoint>();
		if (ppt) val += ppt->m_j.y;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteFluxZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESoluteMaterialPoint* ppt = el.m_State[i]->ExtractData<FESoluteMaterialPoint>();
		if (ppt) val += ppt->m_j.z;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteRefConcentration::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESoluteMaterialPoint* ppt = el.m_State[i]->ExtractData<FESoluteMaterialPoint>();
		if (ppt) val += ppt->m_crc;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteConcentration_::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.m_State[i]->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_ca[m_nsol];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteFluxX_::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.m_State[i]->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_j[m_nsol].x;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteFluxY_::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.m_State[i]->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_j[m_nsol].y;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteFluxZ_::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.m_State[i]->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_j[m_nsol].z;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemElectricPotential::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.m_State[i]->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_psi;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemCurrentDensityX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.m_State[i]->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_Ie.x;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemCurrentDensityY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.m_State[i]->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_Ie.y;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemCurrentDensityZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.m_State[i]->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_Ie.z;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSBMConcentration_::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.m_State[i]->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_sbmr[m_nsol];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
void ElementDataRecord::Parse(const char *szexpr)
{
	char szcopy[MAX_STRING] = {0};
	strcpy(szcopy, szexpr);
	char* sz = szcopy, *ch;
	m_Data.clear();
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
		if      (strcmp(sz, "x"  ) == 0) m_Data.push_back(new FELogElemPosX(m_pfem));
		else if (strcmp(sz, "y"  ) == 0) m_Data.push_back(new FELogElemPosY(m_pfem));
		else if (strcmp(sz, "z"  ) == 0) m_Data.push_back(new FELogElemPosZ(m_pfem));
		else if (strcmp(sz, "J"  ) == 0) m_Data.push_back(new FELogElemJacobian(m_pfem));
		else if (strcmp(sz, "Ex" ) == 0) m_Data.push_back(new FELogElemStrainX(m_pfem));
		else if (strcmp(sz, "Ey" ) == 0) m_Data.push_back(new FELogElemStrainY(m_pfem));
		else if (strcmp(sz, "Ez" ) == 0) m_Data.push_back(new FELogElemStrainZ(m_pfem));
		else if (strcmp(sz, "Exy") == 0) m_Data.push_back(new FELogElemStrainXY(m_pfem));
		else if (strcmp(sz, "Eyz") == 0) m_Data.push_back(new FELogElemStrainYZ(m_pfem));
		else if (strcmp(sz, "Exz") == 0) m_Data.push_back(new FELogElemStrainXZ(m_pfem));
		else if (strcmp(sz, "E1" ) == 0) m_Data.push_back(new FELogElemStrain1(m_pfem));
		else if (strcmp(sz, "E2" ) == 0) m_Data.push_back(new FELogElemStrain2(m_pfem));
		else if (strcmp(sz, "E3" ) == 0) m_Data.push_back(new FELogElemStrain3(m_pfem));
		else if (strcmp(sz, "sx" ) == 0) m_Data.push_back(new FELogElemStressX(m_pfem));
		else if (strcmp(sz, "sy" ) == 0) m_Data.push_back(new FELogElemStressY(m_pfem));
		else if (strcmp(sz, "sz" ) == 0) m_Data.push_back(new FELogElemStressZ(m_pfem));
		else if (strcmp(sz, "sxy") == 0) m_Data.push_back(new FELogElemStressXY(m_pfem));
		else if (strcmp(sz, "syz") == 0) m_Data.push_back(new FELogElemStressYZ(m_pfem));
		else if (strcmp(sz, "sxz") == 0) m_Data.push_back(new FELogElemStressXZ(m_pfem));
		else if (strcmp(sz, "s1" ) == 0) m_Data.push_back(new FELogElemStress1(m_pfem));
		else if (strcmp(sz, "s2" ) == 0) m_Data.push_back(new FELogElemStress2(m_pfem));
		else if (strcmp(sz, "s3" ) == 0) m_Data.push_back(new FELogElemStress3(m_pfem));
		else if (strcmp(sz, "Fxx") == 0) m_Data.push_back(new FELogElemDeformationGradientXX(m_pfem));
		else if (strcmp(sz, "Fyy") == 0) m_Data.push_back(new FELogElemDeformationGradientXY(m_pfem));
		else if (strcmp(sz, "Fzz") == 0) m_Data.push_back(new FELogElemDeformationGradientXZ(m_pfem));
		else if (strcmp(sz, "Fyz") == 0) m_Data.push_back(new FELogElemDeformationGradientYX(m_pfem));
		else if (strcmp(sz, "Fzx") == 0) m_Data.push_back(new FELogElemDeformationGradientYY(m_pfem));
		else if (strcmp(sz, "Fxy") == 0) m_Data.push_back(new FELogElemDeformationGradientYZ(m_pfem));
		else if (strcmp(sz, "Fyx") == 0) m_Data.push_back(new FELogElemDeformationGradientZX(m_pfem));
		else if (strcmp(sz, "Fxz") == 0) m_Data.push_back(new FELogElemDeformationGradientZY(m_pfem));
		else if (strcmp(sz, "Fzy") == 0) m_Data.push_back(new FELogElemDeformationGradientZZ(m_pfem));
		else if (strcmp(sz, "p"  ) == 0) m_Data.push_back(new FELogElemFluidPressure(m_pfem));
		else if (strcmp(sz, "wx" ) == 0) m_Data.push_back(new FELogElemFluidFluxX(m_pfem));
		else if (strcmp(sz, "wy" ) == 0) m_Data.push_back(new FELogElemFluidFluxY(m_pfem));
		else if (strcmp(sz, "wz" ) == 0) m_Data.push_back(new FELogElemFluidFluxZ(m_pfem));
		else if (strcmp(sz, "c"  ) == 0) m_Data.push_back(new FELogElemSoluteConcentration(m_pfem));
		else if (strcmp(sz, "jx" ) == 0) m_Data.push_back(new FELogElemSoluteFluxX(m_pfem));
		else if (strcmp(sz, "jy" ) == 0) m_Data.push_back(new FELogElemSoluteFluxY(m_pfem));
		else if (strcmp(sz, "jz" ) == 0) m_Data.push_back(new FELogElemSoluteFluxZ(m_pfem));
		else if (strcmp(sz, "crc") == 0) m_Data.push_back(new FELogElemSoluteRefConcentration(m_pfem));
		else if (strcmp(sz, "c1"  ) == 0) m_Data.push_back(new FELogElemSoluteConcentration_T<0>(m_pfem));
		else if (strcmp(sz, "c2"  ) == 0) m_Data.push_back(new FELogElemSoluteConcentration_T<1>(m_pfem));
		else if (strcmp(sz, "c3"  ) == 0) m_Data.push_back(new FELogElemSoluteConcentration_T<2>(m_pfem));
		else if (strcmp(sz, "c4"  ) == 0) m_Data.push_back(new FELogElemSoluteConcentration_T<3>(m_pfem));
		else if (strcmp(sz, "c5"  ) == 0) m_Data.push_back(new FELogElemSoluteConcentration_T<4>(m_pfem));
		else if (strcmp(sz, "c6"  ) == 0) m_Data.push_back(new FELogElemSoluteConcentration_T<5>(m_pfem));
		else if (strcmp(sz, "c7"  ) == 0) m_Data.push_back(new FELogElemSoluteConcentration_T<6>(m_pfem));
		else if (strcmp(sz, "c8"  ) == 0) m_Data.push_back(new FELogElemSoluteConcentration_T<7>(m_pfem));
		else if (strcmp(sz, "j1x" ) == 0) m_Data.push_back(new FELogElemSoluteFluxX_T<0>(m_pfem));
		else if (strcmp(sz, "j1y" ) == 0) m_Data.push_back(new FELogElemSoluteFluxY_T<0>(m_pfem));
		else if (strcmp(sz, "j1z" ) == 0) m_Data.push_back(new FELogElemSoluteFluxZ_T<0>(m_pfem));
		else if (strcmp(sz, "j2x" ) == 0) m_Data.push_back(new FELogElemSoluteFluxX_T<1>(m_pfem));
		else if (strcmp(sz, "j2y" ) == 0) m_Data.push_back(new FELogElemSoluteFluxY_T<1>(m_pfem));
		else if (strcmp(sz, "j2z" ) == 0) m_Data.push_back(new FELogElemSoluteFluxZ_T<1>(m_pfem));
		else if (strcmp(sz, "j3x" ) == 0) m_Data.push_back(new FELogElemSoluteFluxX_T<2>(m_pfem));
		else if (strcmp(sz, "j3y" ) == 0) m_Data.push_back(new FELogElemSoluteFluxY_T<2>(m_pfem));
		else if (strcmp(sz, "j3z" ) == 0) m_Data.push_back(new FELogElemSoluteFluxZ_T<2>(m_pfem));
		else if (strcmp(sz, "j4x" ) == 0) m_Data.push_back(new FELogElemSoluteFluxX_T<3>(m_pfem));
		else if (strcmp(sz, "j4y" ) == 0) m_Data.push_back(new FELogElemSoluteFluxY_T<3>(m_pfem));
		else if (strcmp(sz, "j4z" ) == 0) m_Data.push_back(new FELogElemSoluteFluxZ_T<3>(m_pfem));
		else if (strcmp(sz, "j5x" ) == 0) m_Data.push_back(new FELogElemSoluteFluxX_T<4>(m_pfem));
		else if (strcmp(sz, "j5y" ) == 0) m_Data.push_back(new FELogElemSoluteFluxY_T<4>(m_pfem));
		else if (strcmp(sz, "j5z" ) == 0) m_Data.push_back(new FELogElemSoluteFluxZ_T<4>(m_pfem));
		else if (strcmp(sz, "j6x" ) == 0) m_Data.push_back(new FELogElemSoluteFluxX_T<5>(m_pfem));
		else if (strcmp(sz, "j6y" ) == 0) m_Data.push_back(new FELogElemSoluteFluxY_T<5>(m_pfem));
		else if (strcmp(sz, "j6z" ) == 0) m_Data.push_back(new FELogElemSoluteFluxZ_T<5>(m_pfem));
		else if (strcmp(sz, "j7x" ) == 0) m_Data.push_back(new FELogElemSoluteFluxX_T<6>(m_pfem));
		else if (strcmp(sz, "j7y" ) == 0) m_Data.push_back(new FELogElemSoluteFluxY_T<6>(m_pfem));
		else if (strcmp(sz, "j7z" ) == 0) m_Data.push_back(new FELogElemSoluteFluxZ_T<6>(m_pfem));
		else if (strcmp(sz, "j8x" ) == 0) m_Data.push_back(new FELogElemSoluteFluxX_T<7>(m_pfem));
		else if (strcmp(sz, "j8y" ) == 0) m_Data.push_back(new FELogElemSoluteFluxY_T<7>(m_pfem));
		else if (strcmp(sz, "j8z" ) == 0) m_Data.push_back(new FELogElemSoluteFluxZ_T<7>(m_pfem));
		else if (strcmp(sz, "psi" ) == 0) m_Data.push_back(new FELogElemElectricPotential(m_pfem));
		else if (strcmp(sz, "Iex" ) == 0) m_Data.push_back(new FELogElemCurrentDensityX(m_pfem));
		else if (strcmp(sz, "Iey" ) == 0) m_Data.push_back(new FELogElemCurrentDensityY(m_pfem));
		else if (strcmp(sz, "Iez" ) == 0) m_Data.push_back(new FELogElemCurrentDensityZ(m_pfem));
		else if (strcmp(sz, "sbm1") == 0) m_Data.push_back(new FELogElemSBMConcentration_T<0>(m_pfem));
		else if (strcmp(sz, "sbm2") == 0) m_Data.push_back(new FELogElemSBMConcentration_T<1>(m_pfem));
		else if (strcmp(sz, "sbm3") == 0) m_Data.push_back(new FELogElemSBMConcentration_T<2>(m_pfem));
		else if (strcmp(sz, "sbm4") == 0) m_Data.push_back(new FELogElemSBMConcentration_T<3>(m_pfem));
		else if (strcmp(sz, "sbm5") == 0) m_Data.push_back(new FELogElemSBMConcentration_T<4>(m_pfem));
		else if (strcmp(sz, "sbm6") == 0) m_Data.push_back(new FELogElemSBMConcentration_T<5>(m_pfem));
		else if (strcmp(sz, "sbm7") == 0) m_Data.push_back(new FELogElemSBMConcentration_T<6>(m_pfem));
		else if (strcmp(sz, "sbm8") == 0) m_Data.push_back(new FELogElemSBMConcentration_T<7>(m_pfem));
		else throw UnknownDataField(sz);
		sz = ch;
	}
	while (ch);
}

//-----------------------------------------------------------------------------
double ElementDataRecord::Evaluate(int item, int ndata)
{
	// make sure we have an ELT
	if (m_ELT.empty()) BuildELT();

	// find the element
	FEMesh& mesh = m_pfem->GetMesh();
	assert((item >= 1) && (item <= mesh.Elements()));
	ELEMREF e = m_ELT[item-1];
	assert((e.ndom != -1) && (e.nid != -1));
	FEElement* pe = &mesh.Domain(e.ndom).ElementRef(e.nid); assert(pe);

	// get the element value
	return m_Data[ndata]->value(*pe);
}

//-----------------------------------------------------------------------------
void ElementDataRecord::BuildELT()
{
	int i, j;
	m_ELT.clear();
	FEMesh& m = m_pfem->GetMesh();
	int NE = m.Elements();
	m_ELT.resize(NE);
	for (i=0; i<NE; ++i) 
	{
		m_ELT[i].ndom = -1;
		m_ELT[i].nid  = -1;
	}

	for (i=0; i<m.Domains(); ++i)
	{
		FEDomain& d = m.Domain(i);
		int ne = d.Elements();
		for (j=0; j<ne; ++j)
		{
			FEElement& el = d.ElementRef(j);
			m_ELT[el.m_nID-1].ndom = i;
			m_ELT[el.m_nID-1].nid  = j;
		}
	}
}

//-----------------------------------------------------------------------------
void ElementDataRecord::SelectAllItems()
{
	FEMesh& m = m_pfem->GetMesh();
	int n = m.Elements();
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = i+1;
}
