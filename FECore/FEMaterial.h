// FEMaterial.h: interface for the FEMaterial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
#define AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "tens4d.h"
#include "LoadCurve.h"		//---> can we delete this?
#include "FEMaterialFactory.h"
#include "FEParameterList.h"
#include "FEMaterialPoint.h"
#include "DumpFile.h"
#include <string.h>

#define INRANGE(x, a, b) ((x)>=(a) && (x)<=(b))
#define IN_RIGHT_OPEN_RANGE(x, a, b) ((x)>=(a) && (x)<(b))

//-----------------------------------------------------------------------------
//! exception to throw during the material initialization phase

class MaterialError
{
public:
	MaterialError(const char* sz, ...);

	const char* Error() { return m_szerr; }

protected:
	char	m_szerr[512];
};

//-----------------------------------------------------------------------------
// Forward declaration of the FEElasticMaterial class. 
// TODO: The only reason I had to do this is to define the FEMaterial::GetElasticMaterial.
// However, this is only a temporary construct so make sure to delete this forward declaration
// when no longer needed.
class FEElasticMaterial;

//-----------------------------------------------------------------------------
//! exception to throw during material initialization phase
class MaterialRangeError
{
public:
	// szvar = name of variable
	// vmin  = inf value
	// vmax  = sup value
	// bl    = inf is allowed
	// br    = sup is allowed
	MaterialRangeError(const char* szvar, double vmin, double vmax, bool bl, bool br) : m_szvar(szvar), m_vmin(vmin), m_vmax(vmax), m_bl(bl), m_br(br) {}

public:
	const char*	m_szvar;
	double	m_vmin, m_vmax;
	bool	m_bl, m_br;
};

//-----------------------------------------------------------------------------
//! Abstract base class for material types

//! From this class all other material classes are derived.

class FEMaterial : public FEParamContainer
{
public:
	FEMaterial();
	virtual ~FEMaterial();

	//! set/get material name
	void SetName(const char* sz) { strcpy(m_szname, sz); }
	const char* GetName() { return m_szname; }

	//! returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() { return 0; };

	//! performs initialization and parameter checking
	virtual void Init(){}

	int GetID() { return m_nID; }
	void SetID(int nid) { m_nID = nid; }

	//! Serialize material data to archive
	virtual void Serialize(DumpFile& ar);

	//! TODO: I need to move this function up the hierarchy 
	//        once I redesign the material library
	virtual FEElasticMaterial* GetElasticMaterial() { return 0; }

protected:
	char	m_szname[128];	//!< name of material
	int		m_nID;			//!< material ID
};

//-----------------------------------------------------------------------------
//! Base class for solid-materials.
//! These materials need to define the stress and tangent functions.
//!
class FESolidMaterial : public FEMaterial
{
public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) = 0;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) = 0;

	//! calculate derivative of stress w.r.t. solute concentration at material point
	mat3ds Tangent_Concentration(FEMaterialPoint& pt, const int isol);
	
	//! return the material density
	virtual double Density() = 0;

	//! return the molar mass
	virtual double MolarMass() = 0;
};

//-----------------------------------------------------------------------------
//! Base class for nested materials. A nested material describes whose material 
//! response depends on the formulation of another user specified material. For
//! instance, the FEViscoElastic is an example of nested materials.
//!
//! NOTE:	The visco-elastic material uses a new format that no longer requires
//!			the m_nBaseMat parameter anymore. In the interest of backward
//!			compatibility I am implementing a hack that allows both formulations.
//!			If the m_nBaseMat is -1, then the new format
//!			is used where the m_pBase is owned by the FENestedMaterial. If
//!			it is not -1, then it uses the old format where m_nBaseMat
//!			points to another material.

class FENestedMaterial : public FESolidMaterial
{
public:
	FENestedMaterial() { m_nBaseMat = -1; m_pBase = 0; }
	virtual ~FENestedMaterial(){}

	//! return solid component's density
	double Density () { return m_pBase->Density(); }

	//! return solid component's molar mass
	double MolarMass () { return m_pBase->MolarMass(); }

	//! serialization
	void Serialize(DumpFile& ar);

	//! Get the elastic component
	FEElasticMaterial* GetElasticMaterial() { return m_pBase->GetElasticMaterial(); }

public:
	int					m_nBaseMat;	//!< material ID of base material (one-based!)
	FESolidMaterial*	m_pBase;	//!< pointer to base material

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Base class for materials that define multiple sub-components, such as
//! biphasic, biphasic-solute, triphasic and solute material classes.
class FEMultiMaterial : public FEMaterial
{
	class Property
	{
	public:
		Property(const char* szname, int nid) { m_szname = szname; m_nid = nid; }
		virtual ~Property() {}

	public:
		int GetID() { return m_nid; }
		const char* GetName() { return m_szname; }

	public:
		virtual FEMaterial* GetMaterial() = 0;
		virtual bool SetMaterial(FEMaterial* pm) = 0;

	protected:
		int				m_nid;
		const char*		m_szname;
	};

	template <class T> class Property_T : public Property
	{
	public:
		Property_T(T** ppm, const char* szname, int nid) : Property(szname, nid), m_ppm(ppm) { *ppm = 0; }
		bool SetMaterial(FEMaterial* pm) { *m_ppm = dynamic_cast<T*>(pm); return ((*m_ppm) != 0); }
		FEMaterial* GetMaterial() { return *m_ppm; }
	private:
		T**	m_ppm;
	};

public:
	FEMultiMaterial(){}
	~FEMultiMaterial(){}

	// return nr of material properties
	int Components() { return (int) m_Mat.size(); }

	// add component (use this in derived material's constructor)
	template<class T> void AddComponent(T** ppm, const char* szname, int nid = 0) { m_Mat.push_back(new Property_T<T>(ppm, szname, nid)); }

	// get/set component attributes
	int FindComponent(const char* sz, int nid = 0);
	FEMaterial* GetComponent(int n) { return m_Mat[n]->GetMaterial(); }
	const char* GetComponentName(int n) { return m_Mat[n]->GetName(); }
	bool SetComponent(int n, FEMaterial* pm) { return m_Mat[n]->SetMaterial(pm); }

protected:
	vector<Property*>	m_Mat;
};

//-----------------------------------------------------------------------------
//! Global solute data
//! This structure uniquely identifies a solute in multiphasic problems
// TODO: Move this to a different file
struct FESoluteData {
	int					m_nID;			//!< solute ID
	char				m_szname[128];	//!< solute name
	double				m_z;			//!< solute charge number
};


#endif // !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
