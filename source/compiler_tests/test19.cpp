#include "cddefines.h"
#include "container_classes.h"

void test(multi_arr<long,3>& arr)
{
	ml3ci p = arr.ptr(0,0,3);
	(p-2)[2] = 1;
}
