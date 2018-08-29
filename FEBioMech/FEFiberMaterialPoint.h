#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Material point data for elastic fiber materials
class FEFiberMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEFiberMaterialPoint(FEMaterialPoint *pt = 0);
    
	//! copy material point data
	FEMaterialPoint* Copy();
    
	//! Serialize data to archive
	void Serialize(DumpStream& ar);
    
public:
    vec3d   m_n0;    //!< fiber direction, reference configuration, global CS
};
