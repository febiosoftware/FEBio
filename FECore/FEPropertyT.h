#pragma once

//-----------------------------------------------------------------------------
template <class T>
class FEPropertyT : public FEProperty
{
private:
	T**	m_pc;	//!< pointer to pointer of property

public:
	FEPropertyT(T** ppc) : FEProperty(T::classID()) { m_pc = ppc; }

	bool IsArray() const override { return false; }
	bool IsType(FECoreBase* pc) const override { return (dynamic_cast<T*>(pc) != nullptr); }
	void SetProperty(FECoreBase* pc) override 
	{ 
		*m_pc = dynamic_cast<T*>(pc);
		pc->SetParent(GetParent());
	}
	int size() const override { return ((*m_pc) == 0 ? 0 : 1); }

	FECoreBase* get(int i) override { return *m_pc; }
	FECoreBase* get(const char* szname) override
	{
		if ((*m_pc)->GetName() == std::string(szname))
			return *m_pc;
		else
			return 0;
	}

	FECoreBase* getFromID(int nid) override
	{
		if (m_pc && (*m_pc) && ((*m_pc)->GetID() == nid)) return *m_pc; else return 0;
	}

	void Serialize(DumpStream& ar) override
	{
		// TODO: Implement this
	}

	bool Init() override
	{
		if (m_pc && (*m_pc)) { return (*m_pc)->Init(); }
		return (IsRequired() == false);
	}

	bool Validate() override
	{
		if (m_pc && (*m_pc)) return (*m_pc)->Validate();
		return true;
	}
};

//-----------------------------------------------------------------------------
//! Use this class to define array material properties
template<class T> class FEVecPropertyT : public FEProperty
{
private:
	typedef std::vector<T*>	Y;
	Y*	m_pmp;		//!< pointer to actual material property

public:
	FEVecPropertyT(Y* p) : FEProperty(T::classID()) { m_pmp = p; }
	T* operator [] (int i) { return (*m_pmp)[i]; }
	const T* operator [] (int i) const { return (*m_pmp)[i]; }

	virtual bool IsArray() const { return true; }
	virtual bool IsType(FECoreBase* pc) const { return (dynamic_cast<T*>(pc) != 0); }
	virtual void SetProperty(FECoreBase* pc) {
		m_pmp->push_back(dynamic_cast<T*>(pc)); 
		pc->SetParent(GetParent());
	}
	virtual int size() const { return (int)m_pmp->size(); }
	virtual FECoreBase* get(int i) { return (*m_pmp)[i]; }

	virtual FECoreBase* get(const char* szname)
	{ 
		std::string name(szname);
		for (int i=0; i<(int) m_pmp->size(); ++i)
		{
			T* p = (*m_pmp)[i];
			if (p->GetName() == name) return p;
		}
		return 0;
	}

	virtual FECoreBase* getFromID(int nid)
	{
		for (int i = 0; i<(int)m_pmp->size(); ++i)
		{
			T* p = (*m_pmp)[i];
			if (p && (p->GetID() == nid)) return p;
		}
		return 0;
	}

	void AddProperty(FECoreBase* pc) { 
		m_pmp->push_back(dynamic_cast<T*>(pc)); 
		pc->SetParent(GetParent());
	}

	void Clear()
	{
		for (int i=0; i<(int) m_pmp->size(); ++i) delete (*m_pmp)[i];
		m_pmp->clear();
	}

	void Insert(int n, T* pc)
	{
		m_pmp->insert(m_pmp->begin()+n, pc);
		pc->SetParent(GetParent());
	}

	void Serialize(DumpStream& ar)
	{
		if (ar.IsSaving())
		{
			int n = size();
			ar << n;
			for (int i = 0; i<n; ++i)
			{
				T* pm = (*m_pmp)[i];
				Write(ar, pm);
			}
		}
		else
		{
			int n = 0;
			ar >> n;
			m_pmp->assign(n, nullptr);
			for (int i = 0; i<n; ++i)
			{
				(*m_pmp)[i] = dynamic_cast<T*>(Read(ar));
			}
		}
	}

	bool Init() {
		if (m_pmp->empty() && IsRequired()) return false;
		for (size_t i = 0; i<m_pmp->size(); ++i)
		{
			if ((*m_pmp)[i])
			{
				if ((*m_pmp)[i]->Init() == false) return false;
			}
			else return false;
		}
		return true;
	}

	bool Validate() {
		if (m_pmp->empty()) return true;
		for (size_t i = 0; i<m_pmp->size(); ++i)
		{
			if ((*m_pmp)[i])
			{
				if ((*m_pmp)[i]->Validate() == false) return false;
			}
		}
		return true;
	}
};

