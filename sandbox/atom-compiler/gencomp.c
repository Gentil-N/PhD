#include "gencomp.h"
#include "cutils.h"

static inline utype get_atom_affiliated_group(struct Pipeline *pipeline, utype stage_id, utype atom_id)
{
    utype arrp_affiliated_group = (pipeline->stage_count * 3 + 1) + pipeline->data[stage_id] + pipeline->data[pipeline->stage_count + stage_id] + atom_id;
    return pipeline->data[arrp_affiliated_group];
}

static inline struct vec2u *get_atom_relative_pos(struct Pipeline *pipeline, utype gene_id, utype stage_id, utype atom_id)
{
    utype arrp_atom_relative_pos = pipeline->global_data_size + pipeline->unique_gene_size * gene_id + pipeline->grid_config.ria_count.x * pipeline->grid_config.ria_count.y + 2 * stage_id * pipeline->atom_count + 2 * atom_id;
    return (struct vec2u*)&pipeline->data[arrp_atom_relative_pos];
}

static inline struct vec2u *get_group_relative_pos(struct Pipeline *pipeline, utype gene_id, utype stage_id, utype atom_id)
{
    utype arrp_group_relative_pos = pipeline->global_data_size + pipeline->unique_gene_size * gene_id + pipeline->grid_config.ria_count.x * pipeline->grid_config.ria_count.y + 2 * pipeline->stage_count * pipeline->atom_count + pipeline->data[pipeline->stage_count * 3] * pipeline->grid_config.ria_size.x * pipeline->grid_config.ria_size.y + pipeline->data[pipeline->stage_count * 2 + stage_id] * 2 + get_atom_affiliated_group(pipeline, stage_id, atom_id) * 2;
    return (struct vec2u*)&pipeline->data[arrp_group_relative_pos];
}

static inline struct vec2u get_atom_absolute_pos(struct Pipeline *pipeline, utype gene_id, utype stage_id, utype atom_id)
{
    struct vec2u *atom_relative_pos = get_atom_relative_pos(pipeline, gene_id, stage_id, atom_id);
    struct vec2u *group_relative_pos = get_group_relative_pos(pipeline, gene_id, stage_id, atom_id);
    return (struct vec2u){atom_relative_pos->x + group_relative_pos->x, atom_relative_pos->y + group_relative_pos->y};
}

void create_pipeline(struct Pipeline *pipeline, const struct GridConfig *grid_config, utype atom_count, utype max_gene, utype max_mutation_per_gene, utype max_atom_mutation, utype max_group_mutation, utype stage_count, const struct StageInfo *stages)
{
    pipeline->grid_config = *grid_config;
    pipeline->atom_count = atom_count;
    pipeline->stage_count = stage_count;
    pipeline->max_gene = max_gene;
    pipeline->max_mutation_per_gene = max_mutation_per_gene;
    pipeline->max_atom_mutation = max_atom_mutation;
    pipeline->max_group_mutation = max_group_mutation;
    pipeline->global_data_size = 0;
    pipeline->global_data_size += sizeof(utype) * (stage_count * 3 + 1); // stage offsets + stage group counts + stage group count offsets + total group count
    utype stage_sizes[stage_count]; // COPY THOSE DATA!!!
    for_loop(i, stage_count)
    {
        stage_sizes[i] = 0;
        stage_sizes[i] += sizeof(utype) * stages[i].group_count; // atom_count_per_group
        stage_sizes[i] += sizeof(utype) * atom_count * 2; // atom_affiliations + group data
        pipeline->global_data_size += stage_sizes[i];

        /// check sum
        utype atom_count_check = 0;
        for_loop(j, stages[i].group_count)
        {
            atom_count_check += stages[i].atom_count_per_group[j];
        }
        assert(atom_count_check == atom_count);
    }
    pipeline->unique_gene_size = 0;
    pipeline->unique_gene_size += sizeof(utype) * grid_config->ria_count.x * grid_config->ria_count.y * stage_count; // ria grid
    pipeline->unique_gene_size += sizeof(utype) * atom_count * 2 * stage_count; // atom relative pos
    for_loop(i, stage_count)
    {
        pipeline->unique_gene_size += sizeof(utype) * grid_config->ria_size.x * grid_config->ria_size.y * stages[i].group_count; // atom relative grid
        pipeline->unique_gene_size += sizeof(utype) * 2 * stages[i].group_count; // group relative pos
    }
    pipeline->unique_mutation_size = 0;
    pipeline->unique_mutation_size += sizeof(utype) * max_group_mutation * (1 /* group's ID */ + 2 * stage_count /* xy's list */) * 2 /* in case of swap */;
    pipeline->unique_mutation_size += sizeof(utype) * max_atom_mutation * (1 /* group's ID */ + 2 * stage_count /* xy's list */) * 2 /* in case of swap */;
    pipeline->data = mem_alloc(pipeline->global_data_size + pipeline->unique_gene_size * pipeline->max_gene + pipeline->unique_mutation_size * pipeline->max_mutation_per_gene * pipeline->max_gene);

    /// Copy global data
    for_loop(i, stage_count)
    {

    }
}

void destroy_pipeline(struct Pipeline *pipeline)
{
    mem_free(pipeline->data);
}
