#include "FEPlotHeatFlux.h"
#include "FEHeatSolidDomain.h"
#include "FEHeatTransferMaterial.h"

//-----------------------------------------------------------------------------
bool FEPlotHeatFlux::Save(FEDomain &dom, FEDataStream& a)
{
	FEHeatSolidDomain* pbd = dynamic_cast<FEHeatSolidDomain*>(&dom);
	if (pbd)
	{
		// loop over all elements
		for (int i=0; i<pbd->Elements(); ++i)
		{
			// get the next element
			FESolidElement& el = pbd->Element(i);

			// calculate average heat flux
			vec3d ew = vec3d(0,0,0);
			for (int j=0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FEHeatMaterialPoint* pt = (mp.ExtractData<FEHeatMaterialPoint>());
				if (pt) ew += pt->m_q;
			}
			ew /= el.GaussPoints();

			// store to buffer
			a << ew;
		}
		return true;
	}

	return false;
}
