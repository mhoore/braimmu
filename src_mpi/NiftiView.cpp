#include "NiftiView.h"

#include "nifti1.h"
#include "nifti1_io.h"

#include <utility>
#include <cstdio>

template<class T>
std::unique_ptr<NiftyView<T>> NiftyView<T>::fromNIM(nifti_image *nim) {
	switch (nim->datatype)
	{
		case DT_UINT8:
			return std::make_unique<NiftyViewSlave<T,uint8_t>(nim->data);
		case DT_UINT16:
			return std::make_unique<NiftyViewSlave<T,uint16_t>(nim->data);
		case DT_UINT32:
			return std::make_unique<NiftyViewSlave<T,uint32_t>(nim->data);
		case DT_FLOAT64:
			return std::make_unique<NiftyViewSlave<T,float>(nim->data);
		case DT_FLOAT32:
			return std::make_unique<NiftyViewSlave<T,double>(nim->data);
		default:
			printf("Error: nifti file data type cannot be read. datatype=%i . \n", nim->datatype);
		exit(1);
	}
}
