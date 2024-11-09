#ifndef __CRAND_H__
#define __CRAND_H__

#include <sys/types.h>

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

void crand_init();

void crand_add_range(uint low, uint high);

uint crand_gen(size_t range_id);

uint crand_generate(uint low, uint high);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif //__CRAND_H__
