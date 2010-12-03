#pragma once

#include <list>
using namespace std;

// forward declaration of the FEPlotData class
class FEPlotData;

class FEPlotDataFactory
{
public:
	class ClassDescriptor
	{
	public:
		ClassDescriptor(const char* szname)
		{
			m_szname = szname;
			FEPlotDataFactory::RegisterPlotData(this);
		}
		virtual FEPlotData* Create() = 0;
		const char* GetName() { return m_szname; }
	protected:
		const char*	m_szname;
	};

	template <typename T> class ClassDescriptor_T : public ClassDescriptor
	{
	public:
		ClassDescriptor_T(const char* szname) : ClassDescriptor(szname){}
		FEPlotData* Create() { return new T; }
	};

public:
	static FEPlotDataFactory* GetInstance();
	static void RegisterPlotData(ClassDescriptor* pcd)	{ GetInstance()->m_list.push_back(pcd); }
	static FEPlotData* Create(const char* sz);

private:
	list<ClassDescriptor*>	m_list;
	static	FEPlotDataFactory*	m_pThis;

private:
	FEPlotDataFactory(){}
};

#define REGISTER_PLOTDATA(theClass, theName) \
	FEPlotDataFactory::ClassDescriptor_T<theClass> _CD_##theClass##_(theName);
