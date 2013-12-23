#include "stdafx.h"
#include "FESolidMaterial.h"

FESolidMaterial::FESolidMaterial(FEModel* pfem) : FEMaterial(pfem) {}

//! return the material density
double FESolidMaterial::Density() { return m_density; }
