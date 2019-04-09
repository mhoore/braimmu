#pragma once

#include <memory>

template<class T>
class NiftyView
{
	public:

	virtual ~NiftyView() {}

	virtual T get(size_t i) = 0;

	static std::unique_ptr<NiftyView<T>> fromNIM(class nifti_image *nim);
};

template<class T, class From>
class NiftyViewSlave : public NiftyView
{
	From* m_ptr;

	public:

		NiftyViewSlave(void* ptr) : m_ptr(ptr) {}

	T get(size_t i) override {
		return m_ptr[i];
	}
};
