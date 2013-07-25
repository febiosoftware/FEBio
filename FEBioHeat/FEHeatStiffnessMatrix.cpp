#include "FEHeatStiffnessMatrix.h"
#include <FECore/FEAnalysis.h>
#include <FECore/FEModel.h>
using namespace FECore;

//-----------------------------------------------------------------------------
FEHeatStiffnessMatrix::FEHeatStiffnessMatrix(SparseMatrix* pK) : FEGlobalMatrix(pK)
{
}

//-----------------------------------------------------------------------------
FEHeatStiffnessMatrix::~FEHeatStiffnessMatrix()
{
}

//-----------------------------------------------------------------------------
bool FEHeatStiffnessMatrix::Create(FEModel* pfem, int neq, bool breset)
{
	// keep a pointer to the FEM object
	m_pfem = pfem;
	FEModel& fem = *m_pfem;
	FEAnalysis* pstep = fem.GetCurrentStep();
	FEMesh& mesh = fem.GetMesh();

	// The first time we come here we build the "static" profile.
	// This static profile stores the contribution to the matrix profile
	// of the "elements" that do not change. Most elements are static except
	// for instance contact elements which can change connectivity in between
	// calls to the Create() function. Storing the static profile instead of
	// reconstructing it every time we come here saves us a lot of time.
	static SparseMatrixProfile MP;

	// begin building the profile
	build_begin(neq);
	{
		// The first time we are here we construct the "static"
		// profile. This profile contains the contribution from
		// all static elements. A static element is defined as
		// an element that never changes its connectity. This 
		// static profile is stored in the MP object. Next time
		// we come here we simply copy the MP object in stead
		// of building it from scratch.
		if (breset)
		{
			MP.clear();

			vector<int> elm;

			// Add all elements to the profile
			// Loop over all active domains
			for (int nd=0; nd<pstep->Domains(); ++nd)
			{
				FEDomain& d = *pstep->Domain(nd);

				for (int j=0; j<d.Elements(); ++j)
				{
					FEElement& el = d.ElementRef(j);
					if (!el.IsRigid())
					{
						d.UnpackLM(el, elm);
						build_add(elm);
					}
				}
			}
		}
		else
		{
			// copy the old static profile
			*m_pMP = MP;
		}
	}
	// All done! We can now finish building the profile and create 
	// the actual sparse matrix. This is done in the following function
	build_end();

	return true;
}
