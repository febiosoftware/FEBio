#include "stdafx.h"
#include "FEDataGenerator.h"
#include "FEMesh.h"
#include "log.h"

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
	felog.SetMode(Logfile::LOG_FILE_AND_SCREEN);

	// check input
	if (part == 0) return false;
	if (szvar == 0) return false;

	int index = 0;
	char szbuf[256] = { 0 };
	strcpy(szbuf, szvar);
	char* chl = strchr(szbuf, '[');
	if (chl)
	{
		*chl++ = 0;
		char* chr = strrchr(chl, ']');
		if (chr == 0) return false;
		*chr = 0;
		index = atoi(chl);
		if (index < 0) return false;
	}

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
		for (int j = 0; j<nint; ++j)
		{
			// evaluate the spatial position of this gauss point
			vec3d x = el.Evaluate(r, j);

			// find the parameter
			FEMaterialPoint* pt = el.GetMaterialPoint(j);
			FEParam* p = pt->FindParameter(szbuf);
			if (p && (index < p->dim()) && (p->type() == FE_PARAM_DOUBLE))
			{
				*p->pvalue<double>(index) = value(x);
			}
			else return false;
		}
	}

	return true;
}
