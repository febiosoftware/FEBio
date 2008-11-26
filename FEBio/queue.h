// queue.h: interface for the queue class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_QUEUE_H__328DC589_C09F_47F7_9605_8B2320835BFD__INCLUDED_)
#define AFX_QUEUE_H__328DC589_C09F_47F7_9605_8B2320835BFD__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "assert.h"

template <typename T> class queue  
{
public:
	queue(int n)
	{
		m_pd = new T[n];
		m_nmax = n;
		m_nsize = 0;
		m_npush = 0;
		m_npop = -1;
	}
	virtual ~queue()
	{
		delete [] m_pd;
	}

	void Clear()
	{
		m_nsize = 0;
		m_npush = 0;
		m_npop = -1;
	}

	void Pop(T& item)
	{
		// make sure there is something to pop
		assert(m_npop >= 0);

		// copy the item
		item = m_pd[m_npop];

		// move the pop counter
		++m_npop;
		m_npop = m_npop%m_nmax;

		// decrease size counter
		--m_nsize;
		if (m_nsize == 0) m_npop = -1;
	}

	void Push(T& item)
	{
		// make sure the queue is not full
		assert(m_npop != m_npush);

		// push the item
		m_pd[m_npush] = item;

		// move the push counter one over
		++m_npush;
		m_npush = m_npush%m_nmax;
		if (m_nsize == 0) m_npop = 0;

		// increase the size counter
		++m_nsize;
	}

	int Size() { return m_nsize; }


protected:
	T*	m_pd;

	int	m_nmax;		// max size of queue
	int	m_nsize;	// size of queue

	int	m_npop;		// next item to pop
	int m_npush;	// next item to push
};

#endif // !defined(AFX_QUEUE_H__328DC589_C09F_47F7_9605_8B2320835BFD__INCLUDED_)
