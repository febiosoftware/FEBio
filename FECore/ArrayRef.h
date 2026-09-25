#pragma once
#include <memory>
#include <vector>


namespace fecore {

	template <class T>
	class ArrayRef
	{
	public:
		ArrayRef(T* data, size_t size) : m_data(data), m_size(size) {}

		ArrayRef(std::vector<T>& v) : m_data(v.data()), m_size(v.size()) {}

		size_t size() const { return m_size; }

		T& operator[](size_t i) const { return m_data[i]; }

		T* data() const { return m_data; }

	private:
		T* m_data = nullptr;
		size_t m_size = 0;
	};

	template <class T>
	class ConstArrayRef
	{
	public:
		ConstArrayRef(const T* data, size_t size) : m_data(data), m_size(size) {}

		ConstArrayRef(const std::vector<T>& v) : m_data(v.data()), m_size(v.size()) {}

		ConstArrayRef(const ArrayRef<T>& arr) : m_data(arr.data()), m_size(arr.size()) {}

		size_t size() const { return m_size; }

		const T& operator[](size_t i) const { return m_data[i]; }

		const T* data() const { return m_data; }

	private:
		const T* m_data = nullptr;
		size_t m_size = 0;
	};

	template <typename T>
	void zero(ArrayRef<T> arr) { memset(arr.data(), 0, arr.size() * sizeof(T)); }
} // namespace fecore
