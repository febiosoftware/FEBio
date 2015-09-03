#include "stdafx.h"
#include "FEBiphasicSoluteAnalysis.h"
#include "FEBioMech/FERigidMaterial.h"
#include "FEBioMech/FEContactInterface.h"
#include "FEBioMech/FEUncoupledMaterial.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FEBiphasicSoluteAnalysis::FEBiphasicSoluteAnalysis(FEModel* pfem) : FEAnalysis(pfem, FE_POROSOLUTE) {}
