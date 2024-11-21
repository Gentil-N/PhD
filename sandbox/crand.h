#ifndef __CRAND_H__
#define __CRAND_H__

#include <sys/types.h>

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

void crand_init();

uint crand_generate(uint low, uint high);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif //__CRAND_H__
