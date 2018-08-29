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
