#pragma once
#include <memory>

#if defined(__CUDACC__)
	#define FECORE_CUDA_HOST_DEVICE __host__ __device__
	#define FECORE_CUDA_HOST        __host__
	#define FECORE_CUDA_DEVICE      __device__
#else
	#define FECORE_CUDA_HOST_DEVICE
	#define FECORE_CUDA_HOST
	#define FECORE_CUDA_DEVICE
#endif

namespace fecore {

	template <class T>
	class ArrayRef
	{
	public:
		FECORE_CUDA_HOST_DEVICE
			ArrayRef(T* data, size_t size) : m_data(data), m_size(size) {}

		ArrayRef(std::vector<T>& v) : m_data(v.data()), m_size(v.size()) {}

		FECORE_CUDA_HOST_DEVICE
			size_t size() const { return m_size; }

		FECORE_CUDA_HOST_DEVICE
			T& operator[](size_t i) const { return m_data[i]; }

		FECORE_CUDA_HOST_DEVICE
			T* data() const { return m_data; }

	private:
		T* m_data = nullptr;
		size_t m_size = 0;
	};

	template <class T>
	class ConstArrayRef
	{
	public:
		FECORE_CUDA_HOST_DEVICE
			ConstArrayRef(const T* data, size_t size) : m_data(data), m_size(size) {}

		ConstArrayRef(const std::vector<T>& v) : m_data(v.data()), m_size(v.size()) {}

		FECORE_CUDA_HOST_DEVICE
			size_t size() const { return m_size; }

		FECORE_CUDA_HOST_DEVICE
			const T& operator[](size_t i) const { return m_data[i]; }

		FECORE_CUDA_HOST_DEVICE
			const T* data() const { return m_data; }

	private:
		const T* m_data = nullptr;
		size_t m_size = 0;
	};

	template <typename T>
	void zero(ArrayRef<T> arr) { memset(arr.data(), 0, arr.size() * sizeof(T)); }
} // namespace fecore
