#pragma once
#include "FECore/FEMesh.h"

class FEBoxMesh : public FEMesh  
{
public:
	FEBoxMesh(FEModel* fem);
	virtual ~FEBoxMesh();

	void Create(int nx, int ny, int nz, vec3d r0, vec3d r1, int nhex = FE_HEX8G8);
};
