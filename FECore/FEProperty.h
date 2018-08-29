#pragma once
#include <vector>
#include "DumpStream.h"

//-----------------------------------------------------------------------------
class FEParam;
class FECoreBase;
class ParamString;

//-----------------------------------------------------------------------------
//! First attempt at conceptualizing material properties.
//! A material property is essentially any material that is nested inside
//! another material definition. To faciliate automation of properties, we
//! created an explicit interface for such properties.
//! \todo I'd like to make this available for all FECoreBase classes, not just materials.
class FECORE_API FEProperty
{
public:
	enum Flags
	{
		Required      = 0x01,		// the property is required
		ValueProperty = 0x02		// the property is a "value" property
	};

public:
	//! Name of the property.
	//! Note that the name is not copied so it must point to a static string.
	const char*		m_szname;
	bool			m_brequired;	// true if this flag is required (false if optional). Used in FEMaterial::Init().
	bool			m_bvalue;		// true if this is a "value" property (i.e. a property that has a param of the same name)

public:
	// Set\Get the name of the property
	FEProperty& SetName(const char* sz);
	const char* GetName() const;

public: // these functions have to be implemented by derived classes

	//! helper function for identifying if this is an array property or not
	virtual bool IsArray() const = 0;

	//! see if the pc parameter is of the correct type for this property
	virtual bool IsType(FECoreBase* pc) const = 0;

	//! set the property
	virtual void SetProperty(FECoreBase* pc) = 0;

	//! return the size of the property
	virtual int size() const = 0;

	//! return a specific property by index
	virtual FECoreBase* get(int i) = 0;

	//! return a specific property by name
	virtual FECoreBase* get(const char* szname) = 0;

	//! return a specific property by ID
	virtual FECoreBase* getFromID(int nid) = 0;

	//! serialize property data
	virtual void Serialize(DumpStream& ar) = 0;

	//! initializatoin
	virtual bool Init() = 0;

	//! validation
	virtual bool Validate() = 0;

	//! Get the parent of this property
	FECoreBase* GetParent() { return m_pParent; }

	//! Set the parent of this property
	void SetParent(FECoreBase* parent) { m_pParent = parent; }

protected:
	//! some helper functions for reading, writing properties
	void Write(DumpStream& ar, FECoreBase* pc);
	FECoreBase* Read(DumpStream& ar);

protected:
	// This class should not be created directly
	FEProperty();
	virtual ~FEProperty();

private:
	FECoreBase* m_pParent; //!< pointer to the "parent" material
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

	virtual bool IsArray() const { return false; }
	virtual bool IsType(FECoreBase* pc) const { return (dynamic_cast<T*>(pc) != 0); }
	virtual void SetProperty(FECoreBase* pc) { m_pmp = dynamic_cast<T*>(pc); }
	virtual int size() const { return (m_pmp == 0 ? 0 : 1); }

	virtual FECoreBase* get(int i) { return m_pmp; }
	virtual FECoreBase* get(const char* szname) { if (m_pmp->GetName() == std::string(szname)) return m_pmp; else return 0; }

	virtual FECoreBase* getFromID(int nid)
	{
		if (m_pmp && (m_pmp->GetID() == nid)) return m_pmp; else return 0;	
	}

	void Serialize(DumpStream& ar)
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
		if (m_pmp) { return m_pmp->Init(); }
		return (m_brequired == false);
	}

	bool Validate() {
		if (m_pmp) return m_pmp->Validate();
		return true;
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
	const T* operator [] (int i) const { return m_pmp[i]; }

	virtual bool IsArray() const { return true; }
	virtual bool IsType(FECoreBase* pc) const { return (dynamic_cast<T*>(pc) != 0); }
	virtual void SetProperty(FECoreBase* pc) { m_pmp.push_back(dynamic_cast<T*>(pc)); }
	virtual int size() const { return (int)m_pmp.size(); }
	virtual FECoreBase* get(int i) { return m_pmp[i]; }

	virtual FECoreBase* get(const char* szname)
	{ 
		std::string name(szname);
		for (int i=0; i<(int) m_pmp.size(); ++i)
		{
			T* p = m_pmp[i];
			if (p->GetName() == name) return p;
		}
		return 0;
	}

	virtual FECoreBase* getFromID(int nid)
	{
		for (int i = 0; i<(int)m_pmp.size(); ++i)
		{
			T* p = m_pmp[i];
			if (p && (p->GetID() == nid)) return p;
		}
		return 0;
	}

	void AddProperty(FECoreBase* pc) { m_pmp.push_back(dynamic_cast<T*>(pc)); }

	void Clear()
	{
		for (int i=0; i<(int) m_pmp.size(); ++i) delete m_pmp[i];
		m_pmp.clear();
	}

	void Insert(int n, T* pc)
	{
		m_pmp.insert(m_pmp.begin()+n, pc);
	}

	void Serialize(DumpStream& ar)
	{
		if (ar.IsSaving())
		{
			int n = size();
			ar << n;
			for (int i = 0; i<n; ++i)
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
			for (int i = 0; i<n; ++i)
			{
				m_pmp[i] = dynamic_cast<T*>(Read(ar));
			}
		}
	}

	bool Init() {
		if (m_pmp.empty() && m_brequired) return false;
		for (size_t i = 0; i<m_pmp.size(); ++i)
		{
			if (m_pmp[i])
			{
				if (m_pmp[i]->Init() == false) return false;
			}
			else return false;
		}
		return true;
	}

	bool Validate() {
		if (m_pmp.empty()) return true;
		for (size_t i = 0; i<m_pmp.size(); ++i)
		{
			if (m_pmp[i])
			{
				if (m_pmp[i]->Validate() == false) return false;
			}
		}
		return true;
	}

};
