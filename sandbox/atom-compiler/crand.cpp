#include "crand.h"

#include <random>

static std::mt19937 GENERATOR;
static std::random_device RD;

void crand_init()
{
    GENERATOR.seed(RD()); // mersenne_twister_engine
}


uint crand_generate(uint low, uint high)
{
    return std::uniform_int_distribution<uint>(low, high)(GENERATOR);
}
