/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include <vector>

// template class for storing 2D data
template <typename T> class table
{
public:
	// constructors
	table() : m_rows(0), m_cols(0) {}
	table(int nrows, int ncols) : m_rows(0), m_cols(0) { resize(nrows, ncols); }
	table(const table& t) { m_data = t.m_data; m_rows = t.m_rows; m_cols = t.m_cols; }

	// assignment operator
	table& operator = (const table& t) { m_data = t.m_data; m_rows = t.m_rows; m_cols = t.m_cols; return (*this); }

	// resize table
	void resize(int nrows, int ncols, T def = T(0))
	{
		std::vector<T> tmp(m_data);
		m_data.assign(nrows*ncols, def);

		int nr = (nrows < m_rows ? nrows : m_rows);
		int nc = (ncols < m_cols ? ncols : m_cols);
		for (int i=0; i<nr; ++i)
		{
			for (int j=0; j<nc; ++j) m_data[i*ncols + j] = tmp[i*m_cols + j];
		}

		m_rows = nrows;
		m_cols = ncols;
	}

	// assign a value to the entire table
	void set(const T& v) 
	{ 
		if (m_data.empty() == false) 
			m_data.assign(m_rows*m_cols, v); 
	}

	// get sizes
	int rows() const { return m_rows; }
	int columns() const { return m_cols; }

	// access operator
	const T& operator () (int i, int j) const { return m_data[i*m_cols + j]; }
	T& operator () (int i, int j) { return m_data[i*m_cols + j]; }

private:
	std::vector<T>	m_data;
	int				m_rows, m_cols;
};
