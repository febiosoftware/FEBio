// stack.h: interface for the stack class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_STACK_H__2734BF15_E8B2_4304_A03D_7CE8FC955606__INCLUDED_)
#define AFX_STACK_H__2734BF15_E8B2_4304_A03D_7CE8FC955606__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//-----------------------------------------------------------------------------
//! This class implements a simple stack data structure. 

//! Our stack ckass provides push/pop functionality. 
//! The peek function allows you to query the top of the stack
//! without removing the queried item from the stack.
//! Note that the maximum size of the stack has to be input in the constructor.
//! So make sure you don't push too much.
//! Also, note that the actual stack data is not allocated until the first push.

template <typename T> class stack  
{
public:
	//! constructor
	stack(int nmaxsize)
	{
		assert(nmaxsize > 0);

		m_pd = 0;
		m_nmax  = nmaxsize;
		m_nsize = 0;
	}

	//! virtual destructor
	virtual ~stack()
	{ 
		if (m_pd) delete [] m_pd; 
	}

	//! push an item on the stack
	void Push(T& item)
	{
		if (m_pd == 0) { m_pd = new T[m_nmax]; assert(m_pd); }
		assert(m_nsize < m_nmax);
		m_pd[m_nsize++] = item;
	}

	//! pop an item from the stack
	void Pop(T& item)
	{
		assert(m_nsize > 0);
		item = m_pd[--m_nsize];
	}

	//! peek at the top of the stack
	void Peek(T& item)
	{
		assert(m_nsize > 0);
		item = m_pd[m_nsize-1];
	}

	//! clear the stack
	void Clear()
	{
		m_nsize = 0;
	}

	//! return the current size of the stack
	int Size   () { return m_nsize; }

	//! returns true if the stack is empty
	bool IsEmpty() { return (m_nsize == 0); }

	//! return maximum size of stack
	int MaxSize() { return m_nmax; }

protected:
	T*	m_pd;		//!< data pointer
	int	m_nsize;	//!< current size of stack
	int	m_nmax;		//!< max size of stack
};

#endif // !defined(AFX_STACK_H__2734BF15_E8B2_4304_A03D_7CE8FC955606__INCLUDED_)
