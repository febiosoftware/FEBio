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
#include <stddef.h>

#define INRANGE(x, a, b) ((x)>=(a) && (x)<=(b))
#define IN_RIGHT_OPEN_RANGE(x, a, b) ((x)>=(a) && (x)<(b))

//-----------------------------------------------------------------------------
// forward declaration of some classes
class FEModel;
class FEElement;

//-----------------------------------------------------------------------------
//! helper functions for reporting material errors

bool MaterialError(const char* sz, ...);

//-----------------------------------------------------------------------------
// Forward declaration of the FEElasticMaterial class. 
// TODO: The only reason I had to do this is to define the FEMaterial::GetElasticMaterial.
// However, this is only a temporary construct so make sure to delete this forward declaration
// when no longer needed.
class FEElasticMaterial;

//-----------------------------------------------------------------------------
//! First attempt at conceptualizing material properties.
//! A material property is essentially any material that is nested inside
//! another material definition. To faciliate automation of properties, we
//! created an explicit interface for such properties.
//! \todo I'd like to make this available for all FECoreBase classes, not just materials.
class FEProperty
{
public:
	//! Name of the property.
	//! Note that the name is not copied so it must point to a static string.
	const char*		m_szname;
	bool			m_brequired;	// true if this flag is required (false if optional). Used in FEMaterial::Init().

public:
	// Set\Get the name of the property
	FEProperty& SetName(const char* sz);
	const char* GetName() const;

public: // these functions have to be implemented by derived classes
	
	//! see if the pc parameter is of the correct type for this property
	virtual bool IsType(FECoreBase* pc) const = 0;

	//! set the property
	virtual void SetProperty(FECoreBase* pc) = 0;

	//! return the size of the property
	virtual int size() = 0;

	//! return a specific property
	virtual FECoreBase* get(int i) = 0;

	//! return a parameter of the property
	virtual FEParam* GetParameter(const ParamString& s) = 0;

	//! serialize property data
	virtual void Serialize(DumpFile& ar) = 0;

	//! initializatoin
	virtual bool Init() = 0;

	//! Get the parent of this material (zero if none)
	FEMaterial* GetParent() { return m_pParent; }
    
	//! Set the parent of this material
	void SetParent(FEMaterial* pmat) { m_pParent = pmat; }

protected:
	//! some helper functions for reading, writing properties
	void Write(DumpFile& ar, FEMaterial* pc);
	FEMaterial* Read(DumpFile& ar);

protected:
	// This class should not be created directly
	FEProperty();
	virtual ~FEProperty();

private:
	FEMaterial* m_pParent; //!< pointer to the "parent" material
};

//-----------------------------------------------------------------------------
//! Use this class to acutally define material properties in material classes.
//! Note that the m_pmp member can be zero if the material property is optional.
template<class T> class FEPropertyT : public FEProperty
{
private:
	T*				m_pmp;		//!< pointer to actual material property

public:
	FEPropertyT() { m_pmp = nullptr; }
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
			Write(ar, m_pmp);
		}
		else
		{
			m_pmp = dynamic_cast<T*>(Read(ar));
		}
	}

	bool Init() { 
		if (m_pmp) { m_pmp->Init(); return true;  }
		return (m_brequired==false); 
	}
};

//-----------------------------------------------------------------------------
//! Use this class to define array material properties
template<class T> class FEVecPropertyT : public FEProperty
{
private:
	std::vector<T*>	m_pmp;		//!< pointer to actual material property

public:
	FEVecPropertyT() {}
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
				Write(ar, pm);
			}
		}
		else
		{
			int n = 0;
			ar >> n;
			m_pmp.assign(n, nullptr);
			for (int i=0; i<n; ++i)
			{
				m_pmp[i] = dynamic_cast<T*>(Read(ar));
			}
		}
	}

	bool Init() {
		if (m_pmp.empty() && m_brequired) return false;
		for (size_t i=0; i<m_pmp.size(); ++i)
		{
			if (m_pmp[i]) m_pmp[i]->Init(); else return false;
		}
		return true;
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
	virtual bool Init();

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

protected:
	//! Add a material property
	//! Call this in the constructor of derived classes to 
	//! build the property list
	void AddProperty(FEProperty* pp, const char* sz, bool brequired = true);

public:
	//! Find the index of a material property
	int FindPropertyIndex(const char* sz);

	//! Set a material property
	bool SetProperty(int nid, FECoreBase* pm);

	//! return actual number of properties
	int Properties();

	//! return a material property
	FECoreBase* GetProperty(int i);

	//! return a property (class)
	FEProperty* FindProperty(const char* sz);

	//! return a material parameter
	FEParam* GetParameter(const ParamString& s);

	//! return the number of material properties defined
	int MaterialProperties() { return (int) m_Prop.size(); }

	//! return a material property
	FEProperty* MaterialProperty(int i) { return m_Prop[i]; }

	//! Find a material component by type
	//! This returns the first occurrence of this type
	FEMaterial* FindComponentByType(const char* sztype);

private:
	char	m_szname[128];	//!< name of material
	int		m_nID;			//!< material ID
	int		m_nRB;			//!< rigid body ID (TODO: I hope to remove this sometime)

private:
	FECoordSysMap*	m_pmap;			//!< local material coordinate system
	FEModel*		m_pfem;			//!< pointer to model this material belongs to
	FEMaterial*		m_pParent;		//!< pointer to "parent" material (if any)

	vector<FEProperty*>	m_Prop;		//!< list of material properties
};

#endif // !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
