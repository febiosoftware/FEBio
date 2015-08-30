#include "stdafx.h"
#include "FEMixDomainFactory.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FEMultiphasic.h"
#include "FEBiphasicSolidDomain.h"
#include "FEBiphasicSoluteDomain.h"
#include "FETriphasicDomain.h"
#include "FEMultiphasicDomain.h"

//-----------------------------------------------------------------------------
FEDomain* FEMixDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel* pfem = pmat->GetFEModel();
	FE_Element_Shape eshape = spec.eshape;
	FEDomain* pd = 0;
	if (dynamic_cast<FEBiphasic*>(pmat))
	{
		// biphasic elements
		if ((eshape == ET_HEX8) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_HEX20) || (eshape == ET_HEX27)) pd = new FEBiphasicSolidDomain(pfem);
		else return 0;
	}
	if (dynamic_cast<FEBiphasicSolute*>(pmat))
	{
		// biphasic solute elements
		if ((eshape == ET_HEX8) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_HEX20) || (eshape == ET_HEX27)) pd = new FEBiphasicSoluteDomain(pfem);
		else return 0;
	}
	else if (dynamic_cast<FETriphasic*>(pmat))
	{
		// triphasic elements
		if ((eshape == ET_HEX8) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_HEX20) || (eshape == ET_HEX27)) pd = new FETriphasicDomain(pfem);
		else return 0;
	}
	if (dynamic_cast<FEMultiphasic*>(pmat))
	{
		// multiphasic elements
		if ((eshape == ET_HEX8) || (eshape == ET_PENTA6) || (eshape == ET_TET4) || (eshape == ET_TET10) || (eshape == ET_HEX20) || (eshape == ET_HEX27)) pd = new FEMultiphasicDomain(pfem);
		else return 0;
	}
	if (pd) pd->SetMaterial(pmat);
	return pd;
}
