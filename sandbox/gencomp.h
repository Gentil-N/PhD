#ifndef __GENCOMP_H__
#define __GENCOMP_H__

#define ATOMCOMP_VERSION_MAJOR 0
#define ATOMCOMP_VERSION_MINOR 1

#include <stddef.h>
#include <stdint.h>
#include <sys/types.h>
#include <assert.h>
#include <limits.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

#define utype uint
#define INVALID_UTYPE UINT_MAX

struct vec2i
{
    union
    {
        struct { int x, y; };
        int data[2];
    };
};

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

struct Grid
{
    utype *cells;
};

struct Gene
{
    bool is_alive;
    bool is_survivor;
    utype score;
    struct Grid *ria_grids; // size = stage_count, internal size = grid_config.ria_count.x * grid_config.ria_count.y
    struct vec2u **atom_relative_pos; // size = stage_count, sub size = atom_count
    struct Grid **group_relative_grids; // size = stage_count, sub size = stage_infos[i].group_count, internal size = grid_config.ria_size.x * grid_config.ria_size.y
    struct vec2u **group_pos; // size = stage_count, sub size = stage_infos[i].group_count
};

struct Pipeline
{
    struct GridConfig grid_config;
    utype stage_count;
    utype atom_count;
    utype gene_alive_count;
    utype mutation_count_per_gene;
    utype max_gene;
    struct StageInfo *stage_infos;
    struct Gene *genes; // size = max_gene
    utype *survivor_ids; // size = gene_alive_count
};

struct LogFile;

void create_pipeline(struct Pipeline *pipeline, const struct GridConfig *grid_config, utype atom_count, utype gene_alive_count, utype mutation_count_per_gene, utype stage_count, struct StageInfo *stage_infos);

void destroy_pipeline(struct Pipeline *pipeline);

void clone_gene(struct Pipeline *pipeline, utype original_gene_id, utype target_gene_id);

void mutate_gene(struct Pipeline *pipeline, utype gene_id, utype group_mutation_count_per_stage, utype atom_mutation_count_per_stage, utype max_xshift_group, utype max_xshift_atom, bool clone);

void mutate_all_genes(struct Pipeline *pipeline, utype group_mutation_count_per_stage, utype atom_mutation_count_per_stage, utype max_xshift_group, utype max_xshift_atom, bool clone);

void measure_all_scores(struct Pipeline *pipeline);

utype darwin(struct Pipeline *pipeline, bool show_generation_stats);

void show_gene_to_console(const struct Pipeline *pipeline, utype gene_id, utype stage_id);

void show_gene_state_to_console(const struct Pipeline *pipeline, utype gene_id);

void show_all_gene_states_to_console(const struct Pipeline *pipeline);

struct LogFile *create_log_file(const char *name);

void destroy_log_file(struct LogFile *log_file);

void log_file_write_header(struct LogFile *log_file, const struct Pipeline *pipeline);

void log_file_write_gene(struct LogFile *log_file, const struct Pipeline *pipeline, utype gene_id, const char *selection_rule_name);

utype generate_random_stage(utype atom_count, utype max_atom_per_group, utype *atom_affiliations);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif //__GENCOMP_H__
