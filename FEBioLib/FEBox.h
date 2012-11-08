// FEBox.h: interface for the FEBox class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEBOX_H__1ABC33AE_1143_4836_A943_4AEA8D51704E__INCLUDED_)
#define AFX_FEBOX_H__1ABC33AE_1143_4836_A943_4AEA8D51704E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/FEMesh.h"

class FEBox : public FEMesh  
{
public:
	FEBox();
	virtual ~FEBox();

	void Create(int nx, int ny, int nz, vec3d r0, vec3d r1, int nhex = FE_HEX8G8);
};

#endif // !defined(AFX_FEBOX_H__1ABC33AE_1143_4836_A943_4AEA8D51704E__INCLUDED_)
