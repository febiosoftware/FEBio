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
	// The first FEModel (fem) is the macro-problem, i.e. the model that will generate the callbacks
	// The second FEModel (rve) is the micro-problem that needs to be tracked.
	FERVEProbe(FEModel& fem, FEModel& rve, const char* szfile) : FECallBack(&fem, CB_MAJOR_ITERS | CB_INIT | CB_SOLVED), m_rve(rve), m_file(szfile) 
	{
		m_xplt = 0;
	}

	void Execute(FEModel& fem, int nwhen)
	{
		if (nwhen == CB_INIT)	// initialize the plot file
		{
			// create a plot file
			m_xplt = new FEBioPlotFile(m_rve);
			if (m_xplt->Open(m_rve, m_file.c_str()) == false)
			{
				felog.printf("Failed creating probe.\n\n");
				delete m_xplt; m_xplt = 0;
			}

			// write the initial state
			m_xplt->Write(m_rve);
		}
		else if (nwhen == CB_MAJOR_ITERS)	// store the current state
		{
			// write the deformed state
			if (m_xplt) m_xplt->Write(m_rve);
		}
		else if (nwhen == CB_SOLVED)	// clean up
		{
			if (m_xplt) delete m_xplt;
			m_xplt = 0;
		}
	}

private:
	FEModel&			m_rve;		//!< The RVE model to keep track of
	FEBioPlotFile*		m_xplt;		//!< the actual plot file
	std::string			m_file;		//!< file name
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
				FEMicroMaterialPoint& mmpt = *mp.ExtractData<FEMicroMaterialPoint>();
				FERVEProbe* prve = new FERVEProbe(fem, mmpt.m_rve, p.m_szfile);
			}
			else return fecore_error("Invalid gausspt number for micro-probe %d in material %d (%s)", i+1, m_pMat->GetID(), m_pMat->GetName());
		}
		else return fecore_error("Invalid Element ID for micro probe %d in material %d (%s)", i+1, m_pMat->GetID(), m_pMat->GetName());
	}

	return true;
}
