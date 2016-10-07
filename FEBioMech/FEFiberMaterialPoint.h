#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Material point data for elastic fiber materials
class FEFiberMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEFiberMaterialPoint(FEMaterialPoint *pt = 0) : FEMaterialPoint(pt) {}
    
	//! copy material point data
	FEMaterialPoint* Copy();
    
	//! Initialize material point data
	void Init();
    
	//! Serialize data to archive
	void Serialize(DumpStream& ar);
    
public:
    vec3d   m_n0;    //!< fiber direction, reference configuration, global CS
};
