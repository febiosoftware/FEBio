#pragma once
#include "FESurfaceTraction.h"
#include <FECore/FESurfaceMap.h>
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! The pressure surface is a surface domain that sustains pressure boundary
//! conditions
//!
class FEPressureLoad : public FESurfaceTraction
{
public:
	//! constructor
	FEPressureLoad(FEModel* pfem);

	//! Set the surface to apply the load to
	void SetSurface(FESurface* ps) override;

	//! serialize data
	void Serialize(DumpStream& ar) override;

protected:
	//! calculate traction
	vec3d Traction(const FESurfaceMaterialPoint& mp) override;

	//! calculate stiffness for an element
	void ElementStiffness(FESurfaceElement& el, matrix& ke) override;

	void SymmetricPressureStiffness(FESurfaceElement& el, matrix& ke);
	void UnsymmetricPressureStiffness(FESurfaceElement& el, matrix& ke);

protected:
	FEParamDouble	m_pressure;	//!< pressure value
	bool			m_bsymm;	//!< use symmetric formulation

	DECLARE_FECORE_CLASS();
};
