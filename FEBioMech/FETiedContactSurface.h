#pragma once
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
//! This class describes a contact slave or master surface used for 
//! tied contact

//!	this class is used in contact analyses to describe a contacting
//! surface in a tied contact interface.

class FETiedContactSurface : public FEContactSurface
{
public:
	//! constructor
	FETiedContactSurface(FEModel* pfem) : FEContactSurface(pfem) { m_boffset = false; }

	//! Initializes data structures
	bool Init();

	//! data serialization
	void Serialize(DumpStream& ar);

	//! offset shell surfaces (must be called before Init())
	void SetShellOffset(bool b) { m_boffset = b; }

public:
    void GetContactGap     (int nface, double& pg);
    void GetContactPressure(int nface, double& pg);
    void GetContactTraction(int nface, vec3d& pt);
	void GetNodalContactGap     (int nface, double* gn);
	void GetNodalContactPressure(int nface, double* pn);
	void GetNodalContactTraction(int nface, vec3d* tn);

public:
	vector<vec3d>				m_gap;	//!< gap function at nodes
	vector<FESurfaceElement*>	m_pme;	//!< master element a slave node penetrates
	vector<vec2d>				m_rs;	//!< natural coordinates of slave projection on master element
	vector<vec3d>				m_Lm;	//!< Lagrange multipliers
	vector<vec3d>				m_Tc;	//!< contact forces
	vector<double>				m_off;	//!< offset values (used for shells)

protected:
	bool	m_boffset;		//!< offset shells
};
