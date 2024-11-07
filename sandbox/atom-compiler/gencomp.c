#include "gencomp.h"
#include "cutils.h"

#include <stdint.h>
#include <stdio.h>

static struct vec2u list_index_to_coordinates(utype list_index, utype column_count)
{
    return (struct vec2u){list_index % column_count, (utype)(list_index / column_count)};
}

static utype coordinates_to_list_index(const struct vec2u *coordinates, utype column_count)
{
    return (utype)(coordinates->y * column_count + coordinates->x);
}

static void create_grid(struct Grid *grid, utype cell_count)
{
    grid->cells = mem_malloc(sizeof(utype) * cell_count);
    for_loop(i, cell_count)
    {
        grid->cells[i] = INVALID_UTYPE;
    }
}

static void destroy_grid(struct Grid *grid)
{
    mem_free(grid->cells);
}

static void create_mutation(struct Mutation *mutation, utype stage_count)
{
    mutation->id = INVALID_UTYPE;
    mutation->new_pos = mem_malloc(sizeof(struct vec2u*) * stage_count);
}

static void destroy_mutation(struct Mutation *mutation)
{
    mem_free(mutation->new_pos);
}

static void create_global_mutation(struct GlobalMutation *global_mutation, utype max_atom_mutation, utype max_group_mutation, utype stage_count)
{
    global_mutation->atom_mutation_count = 0;
    global_mutation->atom_mutations = mem_malloc(sizeof(struct Mutation) * max_atom_mutation);
    for_loop(i, max_atom_mutation)
    {
        create_mutation(&global_mutation->atom_mutations[i], stage_count);
    }
    global_mutation->group_mutation_count = 0;
    global_mutation->group_mutations = mem_malloc(sizeof(struct Mutation) * max_group_mutation);
    for_loop(i, max_group_mutation)
    {
        create_mutation(&global_mutation->group_mutations[i], stage_count);
    }
}

static void destroy_global_mutation(struct GlobalMutation *global_mutation, utype max_atom_mutation, utype max_group_mutation)
{
    for_loop(i, max_group_mutation)
    {
        destroy_mutation(&global_mutation->group_mutations[i]);
    }
    mem_free(global_mutation->group_mutations);
    for_loop(i, max_atom_mutation)
    {
        destroy_mutation(&global_mutation->atom_mutations[i]);
    }
    mem_free(global_mutation->atom_mutations);
}

static void create_gene(struct Pipeline *pipeline, struct Gene *gene)
{
    /// allocating resources
    gene->ria_grids = mem_malloc(sizeof(struct Grid) * pipeline->stage_count);
    for_loop(i, pipeline->stage_count)
    {
        create_grid(&gene->ria_grids[i], pipeline->grid_config.ria_count.x * pipeline->grid_config.ria_count.y);
    }
    gene->atom_relative_pos = mem_malloc(sizeof(struct vec2u*) * pipeline->stage_count);
    for_loop(i, pipeline->stage_count)
    {
        gene->atom_relative_pos[i] = mem_malloc(sizeof(struct vec2u) * pipeline->atom_count);
    }
    gene->group_relative_grids = mem_malloc(sizeof(struct Grid*) * pipeline->stage_count);
    for_loop(i, pipeline->stage_count)
    {
        gene->group_relative_grids[i] = mem_malloc(sizeof(struct Grid) * pipeline->stage_infos[i].group_count);
        for_loop(j, pipeline->stage_infos[i].group_count)
        {
            create_grid(&gene->group_relative_grids[i][j], pipeline->grid_config.ria_size.x * pipeline->grid_config.ria_size.y);
        }
    }
    gene->group_pos = mem_malloc(sizeof(struct vec2u*) * pipeline->stage_count);
    for_loop(i, pipeline->stage_count)
    {
        gene->group_pos[i] = mem_malloc(sizeof(struct vec2u) * pipeline->stage_infos[i].group_count);
    }
    gene->global_mutations = mem_malloc(sizeof(struct GlobalMutation) * pipeline->max_mutation_per_gene);
    for_loop(i, pipeline->max_mutation_per_gene)
    {
        create_global_mutation(&gene->global_mutations[i], pipeline->max_atom_mutation, pipeline->max_group_mutation, pipeline->stage_count);
    }
    /// initializing default
    for_loop(i, pipeline->stage_count)
    {
        for_loop(j, pipeline->stage_infos[i].group_count)
        {
            gene->ria_grids[i].cells[j] = (utype)j; // assuming group_count < grid_config.ria_count.x * grid_config.ria_count.y
            gene->group_pos[i][j] = (struct vec2u){ j % pipeline->grid_config.ria_count.x, (utype)(j / pipeline->grid_config.ria_count.x)};
        }
        for_loop(j, pipeline->atom_count)
        {
            utype group_id = pipeline->stage_infos[i].atom_affiliations[j];
            /// find next empty cell for the atom and set its relative position
            for_loop(k, pipeline->grid_config.ria_size.x * pipeline->grid_config.ria_size.y)
            {
                if (gene->group_relative_grids[i][group_id].cells[k] == INVALID_UTYPE)
                {
                    gene->group_relative_grids[i][group_id].cells[k] = j;
                    gene->atom_relative_pos[i][j] = list_index_to_coordinates(k, pipeline->grid_config.ria_size.x);
                    break;
                }
            }
        }
    }
}

static void destroy_gene(struct Pipeline *pipeline, struct Gene *gene)
{
    for_loop(i, pipeline->max_mutation_per_gene)
    {
        destroy_global_mutation(&gene->global_mutations[i], pipeline->max_atom_mutation, pipeline->max_group_mutation);
    }
    mem_free(gene->global_mutations);
    for_loop(i, pipeline->stage_count)
    {
        mem_free(gene->group_pos[i]);
    }
    mem_free(gene->group_pos);
    for_loop(i, pipeline->stage_count)
    {
        for_loop(j, pipeline->stage_infos[i].group_count)
        {
            destroy_grid(&gene->group_relative_grids[i][j]);
        }
        mem_free(gene->group_relative_grids[i]);
    }
    mem_free(gene->group_relative_grids);
    for_loop(i, pipeline->stage_count)
    {
        mem_free(gene->atom_relative_pos[i]);
    }
    mem_free(gene->atom_relative_pos);
    for_loop(i, pipeline->stage_count)
    {
        destroy_grid(&gene->ria_grids[i]);
    }
    mem_free(gene->ria_grids);
}

void create_pipeline(struct Pipeline *pipeline, const struct GridConfig *grid_config, utype atom_count, utype max_gene, utype max_mutation_per_gene, utype max_atom_mutation, utype max_group_mutation, utype stage_count, struct StageInfo *stage_infos)
{
    pipeline->grid_config = *grid_config;
    pipeline->atom_count = atom_count;
    pipeline->stage_count = stage_count;
    pipeline->max_gene = max_gene;
    pipeline->max_mutation_per_gene = max_mutation_per_gene;
    pipeline->max_atom_mutation = max_atom_mutation;
    pipeline->max_group_mutation = max_group_mutation;
    pipeline->stage_infos = stage_infos;
    /// check sum
    for_loop(i, stage_count)
    {
        utype atom_count_check = 0;
        for_loop(j, stage_infos[i].group_count)
        {
            atom_count_check += stage_infos[i].atom_count_per_group[j];
        }
        assert(atom_count_check == atom_count);
    }
    /// init genes
    pipeline->genes = mem_malloc(sizeof(struct Gene) * pipeline->max_gene);
    for_loop(i, pipeline->max_gene)
    {
        create_gene(pipeline, &pipeline->genes[i]);
    }
    printf("Everything has been created!\n");
}

void destroy_pipeline(struct Pipeline *pipeline)
{
    for_loop(i, pipeline->max_gene)
    {
        destroy_gene(pipeline, &pipeline->genes[i]);
    }
    mem_free(pipeline->genes);
    printf("Everything has been destroyed...\n");
}

void show_gene_to_console(struct Pipeline *pipeline, utype stage_id, utype gene_id)
{
    for_loop(j, pipeline->grid_config.ria_count.y * pipeline->grid_config.ria_size.y)
    {
        for_loop(i, pipeline->grid_config.ria_count.x * pipeline->grid_config.ria_size.x)
        {
            struct vec2u ria_coords = { (utype)(i / pipeline->grid_config.ria_size.x), (utype)(j / pipeline->grid_config.ria_size.y)};
            utype ria_list_index = coordinates_to_list_index(&ria_coords, pipeline->grid_config.ria_count.x);
            utype group_id = pipeline->genes[gene_id].ria_grids[stage_id].cells[ria_list_index];
            if (group_id != INVALID_UTYPE)
            {
                struct vec2u relative_coords = { i - ria_coords.x * pipeline->grid_config.ria_size.x, j - ria_coords.y * pipeline->grid_config.ria_size.y };
                utype cell_list_index = coordinates_to_list_index(&relative_coords, pipeline->grid_config.ria_size.x);
                utype atom_id = pipeline->genes[gene_id].group_relative_grids[stage_id][group_id].cells[cell_list_index];
                if (atom_id != INVALID_UTYPE)
                {
                    //printf("i=%u j=%u ID=%u\t", (utype)i, (utype)j, atom_id);
                    printf("%u\t", atom_id);
                }
                else
                {
                    //printf("hi=%u j=%u ID=.\t", (utype)i, (utype)j);
                    printf(".\t");
                }
            }
            else
            {
                //printf("i=%u j=%u ID=.\t", (utype)i, (utype)j);
                printf(".\t");
            }
        }
        printf("\n");
    }
}
