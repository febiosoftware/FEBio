#pragma once

#include "FESurface.h"
#include "vec2d.h"

//-----------------------------------------------------------------------------
//! This class describes a contact slave or master surface

//!	this class is used in contact analyses to describe a contacting
//! surface in a contact interface.

class FEContactSurface : public FESurface
{
public:
	//! constructor
	FEContactSurface(FEMesh* pm=0);

	//! destructor
	~FEContactSurface();

	//! Set the sibling of this contact surface
	void SetSibling(FEContactSurface* ps);

public:
	virtual void GetNodalContactGap     (int nface, double* pg);
	virtual void GetNodalContactPressure(int nface, double* pg);
	virtual void GetNodalContactTraction(int nface, vec3d* pt);

	virtual vec3d GetContactForce();

protected:
	FEContactSurface* m_pSibling;
};
