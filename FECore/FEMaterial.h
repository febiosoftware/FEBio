// FEMaterial.h: interface for the FEMaterial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
#define AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "tens4d.h"
#include "FECoreBase.h"
#include "FEMaterialPoint.h"
#include "FECoordSysMap.h"
#include "DumpFile.h"
#include "FECoreKernel.h"
#include <string.h>

#define INRANGE(x, a, b) ((x)>=(a) && (x)<=(b))
#define IN_RIGHT_OPEN_RANGE(x, a, b) ((x)>=(a) && (x)<(b))

//-----------------------------------------------------------------------------
// forward declaration of some classes
class FEModel;
class FEElement;

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
//! First attempt at conceptualizing material properties
class FEProperty
{
public:
	const char*		m_szname;	//!< name of property
	int				m_nID;		//!< ID of property

public:
	FEProperty& SetName(const char* sz) { m_szname = sz; return *this; }
	FEProperty& SetID(int nid) { m_nID = nid; return *this; }
	const char* GetName() const { return m_szname; }
	const int GetID() const { return m_nID; }

	virtual bool IsType(FECoreBase* pc) const = 0;
	virtual void SetProperty(FECoreBase* pc) = 0;

	virtual void Serialize(DumpFile& ar) = 0;

	virtual int size() = 0;
	virtual FECoreBase* get(int i) = 0;

	virtual FEParam* GetParameter(const ParamString& s) = 0;

protected:
	FEProperty(){}
	virtual ~FEProperty(){}
};

//-----------------------------------------------------------------------------
//! Use this class to acutally define material properties in class
//! Note that the m_pmp member can be zero if the material property is optional
template<class T> class FEPropertyT : public FEProperty
{
private:
	T*				m_pmp;		//!< pointer to actual material property

public:
	FEPropertyT() { m_szname = 0; m_nID = -1; m_pmp = 0; }
	operator T*() { return m_pmp; }
	T* operator->() { return m_pmp; }
	void operator = (T* p) { m_pmp = p; }

	virtual bool IsType(FECoreBase* pc) const { return (dynamic_cast<T*>(pc) != 0); }
	virtual void SetProperty(FECoreBase* pc) { m_pmp = dynamic_cast<T*>(pc); }
	virtual int size() { return (m_pmp == 0 ? 0 : 1); }
	virtual FECoreBase* get(int i) { return m_pmp; }

	FEParam* GetParameter(const ParamString& s)
	{
		return (m_pmp ? m_pmp->GetParameter(s) : 0); 
	}

	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			int nflag = (m_pmp == 0 ? 0 : 1);
			ar << nflag;
			if (nflag)
			{
				ar << m_pmp->GetTypeStr();
				m_pmp->Serialize(ar);
			}
		}
		else
		{
			int nflag = 0;
			ar >> nflag;
			if (nflag)
			{
				char sz[256];
				ar >> sz;
				m_pmp = dynamic_cast<T*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
				m_pmp->Serialize(ar);

				m_pmp->Init();
			}
		}
	}
};

//-----------------------------------------------------------------------------
//! Use this class to define array material properties
template<class T> class FEVecPropertyT : public FEProperty
{
private:
	std::vector<T*>	m_pmp;		//!< pointer to actual material property

public:
	FEVecPropertyT() { m_szname = 0; m_nID = -1; }
	T* operator [] (int i) { return m_pmp[i]; }

	virtual bool IsType(FECoreBase* pc) const { return (dynamic_cast<T*>(pc) != 0); }
	virtual void SetProperty(FECoreBase* pc) { m_pmp.push_back(dynamic_cast<T*>(pc)); }
	virtual int size() { return (int)m_pmp.size(); }
	virtual FECoreBase* get(int i) { return m_pmp[i]; }

	FEParam* GetParameter(const ParamString& s)
	{
		ParamString s2 = s.next();

		int N = (int)m_pmp.size();
		for (int i = 0; i<N; ++i)
		{
			T* pi = m_pmp[i];
			if (s2 == pi->GetName()) return pi->GetParameter(s2.next());
		}
		return 0;
	}

	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			int n = size();
			ar << n;
			for (int i=0; i<n; ++i)
			{
				T* pm = m_pmp[i];
				int nflag = (pm == 0 ? 0 : 1);
				ar << nflag;
				if (nflag) 
				{
					ar << pm->GetTypeStr();
					pm->Serialize(ar);
				}
			}
		}
		else
		{
			int n = 0;
			ar >> n;
			m_pmp.assign(n, nullptr);
			for (int i=0; i<n; ++i)
			{
				int nflag = 0;
				ar >> nflag;
				if (nflag)
				{
					char sz[256];
					ar >> sz;
					m_pmp[i] = dynamic_cast<T*>(fecore_new<FEMaterial>(FEMATERIAL_ID, sz, ar.GetFEModel()));
					m_pmp[i]->Serialize(ar);
				}
			}
		}
	}
};

//-----------------------------------------------------------------------------
//! Abstract base class for material types

//! From this class all other material classes are derived.

class FEMaterial : public FECoreBase
{
public:
	FEMaterial(FEModel* pfem);
	virtual ~FEMaterial();

	//! set material name
	void SetName(const char* sz);

	//! get the material's name
	const char* GetName();

	//! returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() { return 0; };

	//! performs initialization and parameter checking
	virtual void Init();

	int GetID() { return m_nID; }
	void SetID(int nid) { m_nID = nid; }

	//! Serialize material data to archive
	virtual void Serialize(DumpFile& ar);

	//! Return elastic material \todo I need to move this function up the hierarchy once I redesign the material library
	virtual FEElasticMaterial* GetElasticMaterial() { return 0; }

public:
	// TODO: Some rigid body stuff is moved to here to avoid RTTI and the definition of rigid materials in FECore, 
	//       as well as simplify some initialization. I hope someday to refactor this a bit.
	//! is this a rigid material
	virtual bool IsRigid() { return false; }

	//! get the ID of the rigid body this material is assigned to (-1 if not)
	int GetRigidBodyID() { return m_nRB; }

	//! Set the rigid body ID this material is assigned to
	void SetRigidBodyID(int rid) { m_nRB = rid; }

	//! return the density
	//! TODO: This was added here because the rigid bodies need it to determine the COM
	virtual double Density() { return 0.0; }

public:
	//! Set the local coordinate system map
	void SetCoordinateSystemMap(FECoordSysMap* pmap);

	//! Get the local coordinate system
	FECoordSysMap* GetCoordinateSystemMap();

	//! Set the local coordinate for a material point
	virtual void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp);

	//! return the model this material belongs to
	FEModel* GetFEModel();

	//! Get the parent of this material (zero if none)
	FEMaterial* GetParent() { return m_pParent; }
    
    //! Get the ancestor of this material (this if none)
    FEMaterial* GetAncestor() { FEMaterial* mp = GetParent(); return (mp ? mp->GetAncestor() : this); }

	//! Set the parent of this material
	void SetParent(FEMaterial* pmat) { m_pParent = pmat; }

public:
	// NOTE: This is the new interface that materials have to implement

	//! return the number of material property classes
	virtual int MaterialProperties() { return 0; }

	//! This function must be overloaded (for now) by derived materials that define material properties
	virtual FEProperty* GetMaterialProperty(int nid) { return nullptr; }

public:
	//! Find the index of a material property
	int FindPropertyIndex(const char* sz);

	//! Set a material property
	bool SetProperty(int nid, FECoreBase* pm);

	//! return actual number of properties
	int Properties();

	//! return a material property
	FECoreBase* GetProperty(int i);

	//! return a material parameter
	FEParam* GetParameter(const ParamString& s);

private:
	char	m_szname[128];	//!< name of material
	int		m_nID;			//!< material ID
	int		m_nRB;			//!< rigid body ID (TODO: I hope to remove this sometime)

private:
	FECoordSysMap*	m_pmap;			//!< local material coordinate system
	FEModel*		m_pfem;			//!< pointer to model this material belongs to
	FEMaterial*		m_pParent;		//!< pointer to "parent" material (if any)
};

#endif // !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
