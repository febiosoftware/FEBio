#pragma once
#include "FESurface.h"

class FEPressureSurface : public FESurface
{
public:
	FEPressureSurface(FEMesh* pm) : FESurface(pm) {}

	void PressureStiffness(FESurfaceElement& el, matrix& ke);

	//! Calculates external pressure forces
	bool PressureForce(FESurfaceElement& el, vector<double>& fe);

	//! Calculates the linear external pressure forces (ie. non-follower forces)
	bool LinearPressureForce(FESurfaceElement& el, vector<double>& fe);

	void Serialize(FEM& fem, Archive& ar);

};
