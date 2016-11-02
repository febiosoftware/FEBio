#include "stdafx.h"
#include "FEDataGenerator.h"
#include "FEMesh.h"

FEDataGenerator::FEDataGenerator() : FECoreBase(FEDATAGENERATOR_ID)
{
}

FEDataGenerator::~FEDataGenerator()
{
}

bool FEDataGenerator::Init()
{
	return true;
}

bool FEDataGenerator::Apply(FEDomain* part, const char* szvar)
{
	// check input
	if (part == 0) return false;
	if (szvar == 0) return false;

	FEMesh& mesh = *part->GetMesh();

	vec3d r[FEElement::MAX_NODES];
	size_t nsize = part->Elements();
	for (size_t i = 0; i<nsize; ++i)
	{
		FEElement& el = part->ElementRef(i);
		int neln = el.Nodes();
		int nint = el.GaussPoints();

		// get the element's coordinates
		for (int j=0; j<neln; ++j) r[j] = mesh.Node(el.m_node[j]).m_r0;

		// evaluate the Gauss points
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			// evaluate the spatial position of this gauss point
			vec3d x = el.Evaluate(r, j);

			// find the parameter
			FEMaterialPoint* pt = el.GetMaterialPoint(j);
			while (pt)
			{
				FEParameterList& pl = pt->GetParameterList();
				FEParam* p = pl.Find(szvar);
				if (p)
				{
					if ((p->dim() == 1) && (p->type() == FE_PARAM_DOUBLE))
					{
						p->value<double>() = value(x);
					}
					pt = 0;
				}
				else
				{
					pt = pt->Next();
					if (pt == 0) return false;
				}
			}
		}
	}

	return true;
}
