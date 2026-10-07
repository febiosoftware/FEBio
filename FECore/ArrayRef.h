#pragma once
#include <memory>
#include <vector>
#include <assert.h>

namespace fecore {

	enum class MemorySpace { Host, Device };

	template <class T>
	class ArrayRef
	{
	public:
		ArrayRef(T* data, size_t size, MemorySpace space = MemorySpace::Host) : m_data(data), m_size(size), m_space(space) {}

		ArrayRef(std::vector<T>& v) : m_data(v.data()), m_size(v.size()), m_space(MemorySpace::Host) {}

		T* data() const { return m_data; }

		size_t size() const { return m_size; }
		
		MemorySpace space() const { return m_space; }

	public:
		// Access operator (use only on host!)
		T& operator[](size_t i) const { 
			assert(m_space == MemorySpace::Host); 
			return m_data[i]; 
		}

	private:
		T* m_data = nullptr;
		size_t m_size = 0;
		MemorySpace m_space;
	};

	template <class T>
	class ConstArrayRef
	{
	public:
		ConstArrayRef(const T* data, size_t size, MemorySpace space = MemorySpace::Host) : m_data(data), m_size(size), m_space(space) {}

		ConstArrayRef(const std::vector<T>& v) : m_data(v.data()), m_size(v.size()), m_space(MemorySpace::Host) {}

		ConstArrayRef(const ArrayRef<T>& arr) : m_data(arr.data()), m_size(arr.size()), m_space(arr.space()) {}

		const T* data() const { return m_data; }

		size_t size() const { return m_size; }

		MemorySpace space() const { return m_space; }

	public:
		// access operator (use only on host!)
		const T& operator[](size_t i) const { 
			assert(m_space == MemorySpace::Host); 
			return m_data[i]; 
		}

	private:
		const T* m_data = nullptr;
		size_t m_size = 0;
		MemorySpace m_space;
	};

	// utility function to zero out an ArrayRef (only works on host memory)
	template <typename T>
	void zero(ArrayRef<T> arr) { 
		assert(arr.space() == MemorySpace::Host); 
		memset(arr.data(), 0, arr.size() * sizeof(T)); 
	}
} // namespace fecore
