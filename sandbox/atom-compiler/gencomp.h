#ifndef __GENCOMP_H__
#define __GENCOMP_H__

#include <stddef.h>
#include <stdint.h>
#include <sys/types.h>
#include <assert.h>

#include "cutils.h"

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

#define utype uint

LIST(utype);

struct vec2u
{
    union
    {
        struct { utype x, y; };
        utype data[2];
    };
};

struct GridConfig
{
    struct vec2u ria_count;
    struct vec2u ria_size;
};

struct StageInfo
{
    utype group_count;
    utype *atom_count_per_group;
    utype *atom_affiliations;
    utype *groups_data;
};

struct Pipeline
{
    struct GridConfig grid_config;
    utype atom_count;
    utype stage_count;
    size_t global_data_size;
    size_t unique_gene_size;
    size_t unique_mutation_size;
    utype max_gene;
    utype max_mutation_per_gene;
    utype max_atom_mutation;
    utype max_group_mutation;
    utype *data;
};

void create_pipeline(struct Pipeline *pipeline, const struct GridConfig *grid_config, utype atom_count, utype max_gene, utype max_mutation_per_gene, utype max_atom_mutation, utype max_group_mutation, utype stage_count, const struct StageInfo *stages);

void destroy_pipeline(struct Pipeline *pipeline);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif //__GENCOMP_H__
