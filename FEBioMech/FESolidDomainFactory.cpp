#include "stdafx.h"
#include "FESolidDomainFactory.h"
#include "FELinearSolidDomain.h"
#include "FERigidMaterial.h"
#include "FEUncoupledMaterial.h"
#include "FEElasticSolidDomain.h"
#include "FEElasticShellDomain.h"
#include "FEElasticTrussDomain.h"
#include "FERigidSolidDomain.h"
#include "FERigidShellDomain.h"
#include "FERemodelingElasticDomain.h"
#include "FEUDGHexDomain.h"
#include "FEUT4Domain.h"
#include "FELinearElastic.h"
#include "FE3FieldElasticSolidDomain.h"
#include "FEDiscreteSpringDomain.h"
#include "FERemodelingElasticMaterial.h"
#include "FECore/FEDiscreteMaterial.h"
#include "FEMicroMaterial2O.h"
#include "FEElasticDomain2O.h"
#include "FESRIElasticSolidDomain.h"

//-----------------------------------------------------------------------------
FEDomain* FESolidDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel* pfem = pmat->GetFEModel();
	const char* sztype = 0;
	FE_Element_Shape eshape = spec.eshape;
	FE_Element_Type etype = spec.etype;
	if (dynamic_cast<FERigidMaterial*>(pmat))
	{
		// rigid elements
		if ((eshape == ET_HEX8) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_TET15) || (eshape == ET_HEX20) || (eshape == ET_HEX27)) sztype = "rigid-solid";
		else if ((eshape == ET_QUAD4) || (eshape == ET_TRI3)) sztype = "rigid-shell";
		else return 0;
	}
	else if (dynamic_cast<FERemodelingElasticMaterial*>(pmat)) sztype = "remodeling-solid";
	else if (dynamic_cast<FEMicroMaterial2O*          >(pmat)) sztype = "elastic-20-solid";
	else if (dynamic_cast<FESolidMaterial*>(pmat))
	{
		// structural elements
		if (eshape == ET_HEX8)
		{
			// three-field implementation for uncoupled materials
			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field_hex)) sztype = "three-field-solid";
			else
			{
				if (etype == FE_HEX8G1) sztype = "udg-hex";
				else sztype = "elastic-solid";
			}
		}
		else if ((eshape == ET_TET10) || (eshape == ET_TET15))
		{
			if ((etype == FE_TET10G8RI4)||(etype == FE_TET10G4RI1))
			{
				sztype = "sri-solid";
			}
			else
			{
				if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field_tet)) sztype = "three-field-solid";
				else sztype = "elastic-solid";
			}
		}
		else if ((eshape == ET_HEX20) || (eshape == ET_HEX27))
		{
//			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field_hex)) sztype = "three-field-solid";
			sztype = "elastic-solid";
		}
		else if (eshape == ET_TET4)
		{
			if (spec.m_but4) sztype = "ut4-solid";
			else sztype = "elastic-solid";
		}
		else if (eshape == ET_PENTA6) 
		{
			// three-field implementation for uncoupled materials
			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field_hex)) sztype = "three-field-solid";
			else sztype = "elastic-solid";
		}
		else if ((eshape == ET_QUAD4) || (eshape == ET_TRI3))sztype = "elastic-shell";
		else if ((eshape == ET_TRUSS2)) sztype = "elastic-truss";
		else return 0;
	}
	else if (dynamic_cast<FEDiscreteMaterial*>(pmat))
	{
		// TODO: check the element type
		sztype = "discrete-spring";
	}

	if (sztype)
	{
		FEDomain* pd = fecore_new<FEDomain>(FEDOMAIN_ID, sztype, pfem);
		if (pd) pd->SetMaterial(pmat);
		return pd;
	}
	else return 0;
}
