#include "stdafx.h"
#include "FEBioMixData.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FEMultiphasic.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
double FENodeConcentration::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	const int dof_C = m_pfem->GetDOFIndex("concentration", 0);	
	return node.get(dof_C); 
}

//-----------------------------------------------------------------------------
double FELogElemFluidPressure::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEBiphasicMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEBiphasicMaterialPoint>();
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
		FEBiphasicMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEBiphasicMaterialPoint>();
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
		FEBiphasicMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEBiphasicMaterialPoint>();
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
		FEBiphasicMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEBiphasicMaterialPoint>();
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_ca[0];
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_j[0].x;
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_j[0].y;
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_j[0].z;
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_sbmr[0];
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
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
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_sbmr[m_nsol];
	}
	return val / (double) nint;
}


