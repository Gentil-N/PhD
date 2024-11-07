#ifndef __GENCOMP_H__
#define __GENCOMP_H__

#include <stddef.h>
#include <stdint.h>
#include <sys/types.h>
#include <assert.h>
#include <limits.h>

#include "cutils.h"

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

#define utype uint
#define INVALID_UTYPE UINT_MAX

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

struct Mutation
{
    utype id; // group's id or atom's id
    struct vec2u *new_pos; // size = stage_count
};

struct GlobalMutation
{
    utype group_mutation_count;
    struct Mutation *group_mutations; // size = 2 * max_atom_mutation
    utype atom_mutation_count;
    struct Mutation *atom_mutations; // size = 2 * max_group_mutation
};

struct Grid
{
    utype *cells;
};

struct Gene
{
    struct Grid *ria_grids; // size = stage_count, internal size = grid_config.ria_count.x * grid_config.ria_count.y
    struct vec2u **atom_relative_pos; // size = stage_count, sub size = atom_count
    struct Grid **group_relative_grids; // size = stage_count, sub size = stage_infos[i].group_count, internal size = grid_config.ria_size.x * grid_config.ria_size.y
    struct vec2u **group_pos; // size = stage_count, sub size = stage_infos[i].group_count
    struct GlobalMutation *global_mutations; // size = max_mutation_per_gene
};

struct Pipeline
{
    struct GridConfig grid_config;
    utype stage_count;
    utype atom_count;
    utype max_gene;
    utype max_mutation_per_gene;
    utype max_atom_mutation;
    utype max_group_mutation;
    struct StageInfo *stage_infos;
    struct Gene *genes; // size = max_gene
};

void create_pipeline(struct Pipeline *pipeline, const struct GridConfig *grid_config, utype atom_count, utype max_gene, utype max_mutation_per_gene, utype max_atom_mutation, utype max_group_mutation, utype stage_count, struct StageInfo *stage_infos);

void destroy_pipeline(struct Pipeline *pipeline);

void show_gene_to_console(struct Pipeline *pipeline, utype stage_id, utype gene_id);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif //__GENCOMP_H__
