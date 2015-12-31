#include "stdafx.h"
#include "FEElasticMultiscaleDomain1O.h"
#include "FEMicroMaterial.h"
#include "FECore/mat3d.h"
#include "FECore/tens6d.h"
#include <FECore/FECallBack.h>
#include <FEBioPlot/FEBioPlotFile.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
class FERVEProbe : public FECallBack
{
public:
	FERVEProbe(FEElasticMultiscaleDomain1O* pdom) : FECallBack(pdom->GetFEModel(), CB_MAJOR_ITERS | CB_INIT | CB_SOLVED), m_pdom(pdom) 
	{
		m_xplt = 0;
		m_eid = 0;
		m_ngp = 0;
	}

	void Execute(FEModel& fem, int nwhen)
	{
		// get the element
		FESolidElement& el = m_pdom->Element(m_eid);

		// get the material point
		FEMaterialPoint& mp = *el.GetMaterialPoint(m_ngp);
		FEMicroMaterialPoint& mmp = *mp.ExtractData<FEMicroMaterialPoint>();

		// get the rve model
		FEModel& rve = mmp.m_rve;

		if (nwhen == CB_INIT)	// initialize the plot file
		{
			// create a plot file
			m_xplt = new FEBioPlotFile(rve);
			if (m_xplt->Open(rve, "rve.xplt") == false)
			{
				felog.printf("Failed creating probe.\n\n");
				delete m_xplt; m_xplt = 0;
			}

			// write the initial state
			m_xplt->Write(rve);
		}
		else if (nwhen == CB_MAJOR_ITERS)	// store the current state
		{
			// write the deformed state
			if (m_xplt) m_xplt->Write(rve);
		}
		else if (nwhen == CB_SOLVED)	// clean up
		{
			if (m_xplt) delete m_xplt;
			m_xplt = 0;
		}
	}

private:
	FEElasticMultiscaleDomain1O* m_pdom;		//!< domain
	int		m_eid;								//!< element ID (local in domain)
	int		m_ngp;								//!< which gauss-point
	FEBioPlotFile*		m_xplt;					//!< the actual plot file
};

//-----------------------------------------------------------------------------
//! constructor
FEElasticMultiscaleDomain1O::FEElasticMultiscaleDomain1O(FEModel* pfem) : FEElasticSolidDomain(pfem)
{
}

//-----------------------------------------------------------------------------
//! intialize domain
bool FEElasticMultiscaleDomain1O::Initialize(FEModel& fem)
{
	if (FEElasticSolidDomain::Initialize(fem) == false) return false;

	// get the material
	FEMicroMaterial* pmat = dynamic_cast<FEMicroMaterial*>(m_pMat);
	if (m_pMat == 0) return false;

	// get the master RVE
	FEModel& rve = pmat->m_mrve;

	// loop over all elements
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j) 
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			FEMicroMaterialPoint& mmpt = *mp.ExtractData<FEMicroMaterialPoint>();

			// create the material point RVEs
			if (mmpt.m_rve_init == false)
			{
				mmpt.m_rve.CopyFrom(rve);
				mmpt.m_rve.Init();
				mmpt.m_rve_prev.CopyFrom(rve);
				mmpt.m_rve_prev.Init();
				mmpt.m_rve_init = true;
			}
		}
	}

	// create a probe
	FERVEProbe* prve = new FERVEProbe(this);

	return true;
}
