#include "stdafx.h"
#include "FEBodyLoad.h"

//-----------------------------------------------------------------------------
FEBodyLoad::FEBodyLoad(FEModel* pfem) : FECoreBase(FEBODYLOAD_ID), m_pfem(pfem) {}

//-----------------------------------------------------------------------------
FEBodyLoad::~FEBodyLoad(){}

