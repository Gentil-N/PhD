#include "crand.h"

#include <vector>
#include <random>

static std::vector<std::uniform_int_distribution<uint>> RANGES;
static std::mt19937 GENERATOR;
std::random_device RD;

void crand_init()
{
    GENERATOR.seed(RD()); // mersenne_twister_engine seeded with rd()
}

void crand_add_range(uint low, uint high)
{
    RANGES.push_back(std::uniform_int_distribution<uint>(low, high));
}

uint crand_gen(size_t range_id)
{
    return RANGES[range_id](GENERATOR);
}