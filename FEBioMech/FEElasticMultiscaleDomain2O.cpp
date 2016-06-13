#include "stdafx.h"
#include "FEElasticMultiscaleDomain2O.h"
#include "FEMicroMaterial2O.h"
#include "FECore/mat3d.h"
#include "FECore/tens6d.h"
#include <FECore/log.h>


//-----------------------------------------------------------------------------
//! constructor
FEElasticMultiscaleDomain2O::FEElasticMultiscaleDomain2O(FEModel* pfem) : FEElasticSolidDomain2O(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initialize element data
bool FEElasticMultiscaleDomain2O::Initialize()
{
	// initialize base class first
	if (FEElasticSolidDomain2O::Initialize() == false) return false;

	FEModel& fem = *GetFEModel();

	// initialze RVEs
	const int NE = FEElement::MAX_NODES;
	vec3d x0[NE], xt[NE], r0, rt;
	FEMesh& m = *GetMesh();
		
	FEMicroMaterial2O* pmat = dynamic_cast<FEMicroMaterial2O*>(m_pMat);
	FERVEModel2O& rve = pmat->m_mrve;
			
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();
		for (int i=0; i<neln; ++i)
		{
			x0[i] = m.Node(el.m_node[i]).m_r0;
			xt[i] = m.Node(el.m_node[i]).m_rt;
		}

		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) 
		{
			r0 = el.Evaluate(x0, j);
			rt = el.Evaluate(xt, j);

			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			FEMicroMaterialPoint2O& mmpt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();

			// initialize the material point RVE
			// This essentially copies the master RVE to the material point RVE
			mmpt2O.m_rve.Init(rve);
		}
	}

	// initialize surface RVEs
	int nnf = 0;
	int NF = m_surf.Elements();
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& face = m_surf.Element(i);
		int nint = face.GaussPoints();
		for (int n=0; n<nint; ++n, ++nnf)
		{
			for (int k=0; k<2; ++k)
			{
				FEMaterialPoint& mp = *m_surf.GetData(nnf).m_pt[k];
				FEMicroMaterialPoint2O& mmpt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();

				// Initialize the material point RVE
				// This essentially copies the master RVE model to the material points
				mmpt2O.m_rve.Init(rve);
			}
		}
	}

	// create the probes
	int NP = pmat->Probes();
	for (int i=0; i<NP; ++i)
	{
		FEMicroProbe& p = pmat->Probe(i);
		FEElement* pel = FindElementFromID(p.m_neid);
		if (pel)
		{
			int nint = pel->GaussPoints();
			int ngp = p.m_ngp - 1;
			if ((ngp>=0)&&(ngp<nint))
			{
				FEMaterialPoint& mp = *pel->GetMaterialPoint(ngp);
				FEMicroMaterialPoint2O& mmpt = *mp.ExtractData<FEMicroMaterialPoint2O>();
				FERVEProbe* prve = new FERVEProbe(fem, mmpt.m_rve, p.m_szfile);
			}
			else return fecore_error("Invalid gausspt number for micro-probe %d in material %d (%s)", i+1, m_pMat->GetID(), m_pMat->GetName());
		}
		else return fecore_error("Invalid Element ID for micro probe %d in material %d (%s)", i+1, m_pMat->GetID(), m_pMat->GetName());
	}

	return true;
}
