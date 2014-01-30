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

//-----------------------------------------------------------------------------
int FESolidDomainFactory::GetDomainType(const FE_Element_Spec& spec, FEMaterial* pmat)
{
	FE_Element_Shape eshape = spec.eshape;
	FE_Element_Type etype = spec.etype;
	if (dynamic_cast<FERigidMaterial*>(pmat))
	{
		// rigid elements
		if ((eshape == ET_HEX8) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_TET15) || (eshape == ET_HEX20)) return FE_RIGID_SOLID_DOMAIN;
		else if ((eshape == ET_QUAD4) || (eshape == ET_TRI3)) return FE_RIGID_SHELL_DOMAIN;
		else return 0;
	}
/*	else if (dynamic_cast<FELinearElastic*>(pmat))
	{
		if ((eshape == ET_HEX8) || (eshape == ET_PENTA6) || (eshape == ET_TET4)) return FE_LINEAR_SOLID_DOMAIN;
		else return 0;
	}
*/	else if (dynamic_cast<FESolidMaterial*>(pmat))
	{
		// structural elements
		if (eshape == ET_HEX8)
		{
			// three-field implementation for uncoupled materials
			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field)) return FE_3F_SOLID_DOMAIN;
			else
			{
				if (etype == FE_HEX8G1) return FE_UDGHEX_DOMAIN;
				else return FE_SOLID_DOMAIN;
			}
		}
		else if ((eshape == ET_HEX20) || (eshape == ET_TET10) || (eshape == ET_TET15))
		{
			return FE_SOLID_DOMAIN;
		}
		else if (eshape == ET_TET4)
		{
			if (spec.m_but4) return FE_UT4_DOMAIN;
			else return FE_SOLID_DOMAIN;
		}
		else if (eshape == ET_PENTA6) 
		{
			// three-field implementation for uncoupled materials
			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field)) return FE_3F_SOLID_DOMAIN;
			else return FE_SOLID_DOMAIN;
		}
		else if ((eshape == ET_QUAD4) || (eshape == ET_TRI3)) return FE_SHELL_DOMAIN;
		else if ((eshape == ET_TRUSS2)) return FE_TRUSS_DOMAIN;
		else return 0;
	}
	else if (dynamic_cast<FEDiscreteMaterial*>(pmat))
	{
		// TODO: check the element type
		return FE_DISCRETE_DOMAIN;
	}

	else return 0;
}

//-----------------------------------------------------------------------------
FEDomain* FESolidDomainFactory::CreateDomain(int dtype, FEMesh* pm, FEMaterial* pmat)
{
	if (dtype == FE_SOLID_DOMAIN       ) 
	{
		if (dynamic_cast<FERemodelingElasticMaterial*>(pmat)) return new FERemodelingElasticDomain(pm, pmat);
		else return new FEElasticSolidDomain(pm, pmat);
	}
	if (dtype == FE_SHELL_DOMAIN       ) return new FEElasticShellDomain      (pm, pmat);
	if (dtype == FE_TRUSS_DOMAIN       ) return new FEElasticTrussDomain      (pm, pmat);
	if (dtype == FE_RIGID_SOLID_DOMAIN ) return new FERigidSolidDomain        (pm, pmat);
	if (dtype == FE_RIGID_SHELL_DOMAIN ) return new FERigidShellDomain        (pm, pmat);
	if (dtype == FE_UDGHEX_DOMAIN      ) return new FEUDGHexDomain            (pm, pmat);
	if (dtype == FE_UT4_DOMAIN         ) return new FEUT4Domain               (pm, pmat);
	if (dtype == FE_3F_SOLID_DOMAIN    ) return new FE3FieldElasticSolidDomain(pm, pmat);
	if (dtype == FE_LINEAR_SOLID_DOMAIN) return new FELinearSolidDomain       (pm, pmat);
	if (dtype == FE_DISCRETE_DOMAIN    ) return new FEDiscreteSpringDomain    (pm, pmat);
	return 0;
}
