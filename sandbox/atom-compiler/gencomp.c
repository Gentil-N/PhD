#include "gencomp.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cutils.h"
#include "crand.h"

static utype rand_ria_count_x()
{
    return (utype) crand_gen(1);
}

static utype rand_ria_count_y()
{
    return (utype) crand_gen(2);
}

static utype rand_ria_size_x()
{
    return (utype) crand_gen(3);
}

static utype rand_ria_size_y()
{
    return (utype) crand_gen(4);
}

typedef utype(*rand_func)();

static struct vec2u random_2d(utype x_limit, utype y_limit, utype x_except, utype y_except, rand_func rand_x, rand_func rand_y)
{
    utype x_val = rand_x();
    utype y_val = rand_y();
    /*if (x_val == x_except && y_val == y_except)
    {
        if (x_limit == 1 && y_limit > 1) goto random_2d_move_y;
        else if (x_limit > 1 && y_limit == 1) goto random_2d_move_x;
        else if (x_limit > 1 && y_limit > 1)
        {
            utype choice = crand_gen(0);
            if (choice == 0) goto random_2d_move_y;
            else goto random_2d_move_x;
        }
        else goto random_2d_no_move;
random_2d_move_x:
random_2d_move_y:
    }
random_2d_no_move:*/
    if ((x_val == x_except && y_val == y_except) && (x_limit > 1 || y_limit > 1))
    {
        while (x_val == x_except && y_val == y_except)
        {
            //printf("haha\n");
            x_val = rand_x();
            y_val = rand_y();
        }
    }
    return (struct vec2u){ x_val, y_val };
}

static struct vec2u random_2d_ria_count(utype x_limit, utype y_limit, utype x_except, utype y_except)
{
    return random_2d(x_limit, y_limit, x_except, y_except, rand_ria_count_x, rand_ria_count_y);
}

static struct vec2u random_2d_ria_size(utype x_limit, utype y_limit, utype x_except, utype y_except)
{
    return random_2d(x_limit, y_limit, x_except, y_except, rand_ria_size_x, rand_ria_size_y);
}

/// /!\ does not handle the case where x_shift == y_shift == 0
static struct vec2u random_2d_shift(int curr_x, int curr_y, int x_shift, int y_shift, int x_low, int x_up, int y_low, int y_up)
{
    /*assert(x_low <= 0 && x_up >= 0 && y_low <= 0 && y_up >= 0);
    if (x_low == x_up && y_low == y_up) return (struct vec2i){ 0, 0 };
    int x_val = 0, y_val = 0;
    do {
        x_val = crand_generate(x_low, x_up);
        y_val = crand_generate(y_low, y_up);
    } while ((x_val == 0 && y_val == 0) );
    return (struct vec2i){ x_val, y_val };*/
    uint x_val = curr_x, y_val = curr_y;
    do
    {
        x_val = crand_generate(max(x_low, curr_x - x_shift), min(x_up, curr_x + x_shift));
        y_val = crand_generate(max(y_low, curr_y - y_shift), min(y_up, curr_y + y_shift));
    }
    while ((x_val == curr_x && y_val == curr_y) );
    return (struct vec2u){ x_val, y_val };
}

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

static void copy_grid(const struct Grid *target_grid, struct Grid *original_grid, utype cell_count)
{
    memcpy(target_grid->cells, original_grid->cells, sizeof(utype) * cell_count);
}

static void create_gene(struct Pipeline *pipeline, struct Gene *gene, bool alive, bool survivor)
{
    /// allocating resources
    gene->is_alive = alive;
    gene->is_survivor = survivor;
    gene->score = INVALID_UTYPE;
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
    /// initializing default
    if (alive == false) return;
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

void create_pipeline(struct Pipeline *pipeline, const struct GridConfig *grid_config, utype atom_count, utype gene_alive_count, utype mutation_count_per_gene, utype stage_count, struct StageInfo *stage_infos)
{
    pipeline->grid_config = *grid_config;
    pipeline->atom_count = atom_count;
    pipeline->stage_count = stage_count;
    pipeline->gene_alive_count = gene_alive_count;
    pipeline->mutation_count_per_gene = mutation_count_per_gene;
    pipeline->max_gene = pipeline->gene_alive_count * (1 + pipeline->mutation_count_per_gene);
    pipeline->stage_infos = stage_infos;
    /// init random
    crand_init();
    crand_add_range(0, 1); // #0
    crand_add_range(0, pipeline->grid_config.ria_count.x - 1); // #1
    crand_add_range(0, pipeline->grid_config.ria_count.y - 1); // #2
    crand_add_range(0, pipeline->grid_config.ria_size.x - 1); // #3
    crand_add_range(0, pipeline->grid_config.ria_size.y - 1); // #4
    crand_add_range(0, pipeline->atom_count - 1); // #5
    /// check sum
    for_loop(i, stage_count)
    {
        crand_add_range(0, stage_infos[i].group_count - 1); // #6+i
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
        bool flag = (bool)(i % (1 + pipeline->mutation_count_per_gene) == 0);
        create_gene(pipeline, &pipeline->genes[i], flag, flag);
    }
    /// init survirvors
    pipeline->survivor_ids = mem_malloc(sizeof(utype) * pipeline->gene_alive_count);
    /// end!
    printf("Everything has been created!\n");
    fflush(stdout);
}

void destroy_pipeline(struct Pipeline *pipeline)
{
    mem_free(pipeline->survivor_ids);
    for_loop(i, pipeline->max_gene)
    {
        destroy_gene(pipeline, &pipeline->genes[i]);
    }
    mem_free(pipeline->genes);
    printf("Everything has been destroyed...\n");
}

static void clone_gene_from_ptr(struct Pipeline *pipeline, const struct Gene *original_gene, struct Gene *target_gene)
{
    target_gene->is_alive = original_gene->is_alive;
    target_gene->is_survivor = false; // a cloned gene is not a survivor from the previous darwin pass!
    target_gene->score = original_gene->score;
    for_loop(i, pipeline->stage_count)
    {
        copy_grid(&target_gene->ria_grids[i], &original_gene->ria_grids[i], pipeline->grid_config.ria_count.x * pipeline->grid_config.ria_count.y);
        memcpy(target_gene->atom_relative_pos[i], original_gene->atom_relative_pos[i], sizeof(struct vec2u) * pipeline->atom_count);
        for_loop(j, pipeline->stage_infos[i].group_count)
        {
            copy_grid(&target_gene->group_relative_grids[i][j], &original_gene->group_relative_grids[i][j], pipeline->grid_config.ria_size.x * pipeline->grid_config.ria_size.y);
        }
        memcpy(target_gene->group_pos[i], original_gene->group_pos[i], sizeof(struct vec2u) * pipeline->stage_infos[i].group_count);
    }
}

void clone_gene(struct Pipeline *pipeline, utype original_gene_id, utype target_gene_id)
{
    assert(original_gene_id < pipeline->max_gene && target_gene_id < pipeline->max_gene && original_gene_id != target_gene_id);
    clone_gene_from_ptr(pipeline, &pipeline->genes[original_gene_id], &pipeline->genes[target_gene_id]);
}

static void mutate_gene_group(struct Pipeline *pipeline, utype gene_id, utype stage_id, utype max_xshift_group)
{
    assert(gene_id < pipeline->max_gene);
    struct Gene *current_gene = &pipeline->genes[gene_id];

    utype randomly_chosen_group_id = crand_generate(0, pipeline->stage_infos[stage_id].group_count - 1);
    struct vec2u coords_of_random = current_gene->group_pos[stage_id][randomly_chosen_group_id];
    utype ria_list_index_of_random = coordinates_to_list_index(&coords_of_random, pipeline->grid_config.ria_count.x);
    struct vec2u new_pos_target = random_2d_shift((int)coords_of_random.x, (int)coords_of_random.y, 2, 2, 0, (int)pipeline->grid_config.ria_count.x - 1, 0, (int)pipeline->grid_config.ria_count.y - 1);
    utype target_ria_list_index = coordinates_to_list_index(&new_pos_target, pipeline->grid_config.ria_count.x);
    utype target_group_id = current_gene->ria_grids[stage_id].cells[target_ria_list_index];

    current_gene->group_pos[stage_id][randomly_chosen_group_id] = new_pos_target;
    if (target_group_id != INVALID_UTYPE) current_gene->group_pos[stage_id][target_group_id] = coords_of_random;
    current_gene->ria_grids[stage_id].cells[ria_list_index_of_random] = target_group_id;
    current_gene->ria_grids[stage_id].cells[target_ria_list_index] = randomly_chosen_group_id;
}

static void mutate_gene_atom_inside_group(struct Pipeline *pipeline, utype gene_id, utype stage_id, utype max_xshift_atom)
{
    assert(gene_id < pipeline->max_gene);
    struct Gene *current_gene = &pipeline->genes[gene_id];

    utype randomly_chosen_atom_id = crand_generate(0, pipeline->atom_count - 1);
    utype affilated_group_id = pipeline->stage_infos[stage_id].atom_affiliations[randomly_chosen_atom_id];
    struct vec2u coords_of_random = current_gene->atom_relative_pos[stage_id][randomly_chosen_atom_id];
    utype group_list_index_of_random = coordinates_to_list_index(&coords_of_random, pipeline->grid_config.ria_size.x);
    //printf("%i %i %i %i\n", x_low, x_up, y_low, y_up);
    struct vec2u new_pos_target = random_2d_shift((int)coords_of_random.x, (int)coords_of_random.y, max_xshift_atom, max_xshift_atom, 0, (int)pipeline->grid_config.ria_size.x - 1, 0, (int)pipeline->grid_config.ria_size.y - 1);
    utype target_atom_list_index = coordinates_to_list_index(&new_pos_target, pipeline->grid_config.ria_size.x);
    utype target_atom_id = current_gene->group_relative_grids[stage_id][affilated_group_id].cells[target_atom_list_index];

    current_gene->atom_relative_pos[stage_id][randomly_chosen_atom_id] = new_pos_target;
    if (target_atom_id != INVALID_UTYPE) current_gene->atom_relative_pos[stage_id][target_atom_id] = coords_of_random;
    current_gene->group_relative_grids[stage_id][affilated_group_id].cells[group_list_index_of_random] = target_atom_id;
    current_gene->group_relative_grids[stage_id][affilated_group_id].cells[target_atom_list_index] = randomly_chosen_atom_id;
}

static utype get_closest_dead_gene(const struct Pipeline *pipeline, utype gene_id)
{
    //utype bound = max(pipeline->max_gene - gene_id, gene_id);
    for_loop(i, pipeline->max_gene)
    {
        utype id_left = 0;
        if (i + 1 < gene_id) id_left = gene_id - i - 1;
        utype id_right = pipeline->max_gene - 1;
        if (gene_id + i + 1 < pipeline->max_gene - 1) id_right = gene_id + i + 1;
        if (pipeline->genes[id_left].is_alive == false) return id_left;
        if (pipeline->genes[id_right].is_alive == false) return id_right;
    }
    printf("unable to find a dead gene");
    exit(1);
}

/*void clone_and_mutate_gene(struct Pipeline *pipeline, utype gene_id)
{
    assert(gene_id < pipeline->max_gene);
    printf("closest dead gene %u\n", get_closest_dead_gene(pipeline, 0));
    mutate_gene_atom_inside_group(pipeline, gene_id, 0);
}*/

void mutate_gene(struct Pipeline *pipeline, utype gene_id, utype group_mutation_count_per_stage, utype atom_mutation_count_per_stage, utype max_xshift_group, utype max_xshift_atom, bool clone)
{
    assert(gene_id < pipeline->max_gene && pipeline->genes[gene_id].is_alive == true);
    utype actual_gene_id = gene_id;
    //printf("before %u\n", actual_gene_id);
    if (clone == true)
    {
        actual_gene_id = get_closest_dead_gene(pipeline, gene_id);
        //printf("clone to %u\n", actual_gene_id);
        clone_gene(pipeline, gene_id, actual_gene_id);
        pipeline->genes[actual_gene_id].score = INVALID_UTYPE; // reset score to flag and avoid duplicated score measurement
        //printf("clone to %u\n", actual_gene_id);
    }
    //printf("mutating %u\n", actual_gene_id);
    for_loop(i, pipeline->stage_count)
    {
        for_loop(j, group_mutation_count_per_stage)
        {
            mutate_gene_group(pipeline, actual_gene_id, i, max_xshift_group);
        }
        for_loop(j, atom_mutation_count_per_stage)
        {
            mutate_gene_atom_inside_group(pipeline, actual_gene_id, i, max_xshift_atom);
        }
    }
}

void mutate_all_genes(struct Pipeline *pipeline, utype group_mutation_count_per_stage, utype atom_mutation_count_per_stage, utype max_xshift_group, utype max_xshift_atom, bool clone)
{
    for_loop(i, pipeline->max_gene)
    {
        if (pipeline->genes[i].is_survivor == false || pipeline->genes[i].is_alive == false) continue; // do not mutate early mutated genes and dead genes
        //printf("gene %lu\n", i);
        for_loop(j, pipeline->mutation_count_per_gene)
        {
            mutate_gene(pipeline, i, group_mutation_count_per_stage, atom_mutation_count_per_stage, max_xshift_group, max_xshift_atom, clone);
            if (clone == false) break; // only one pass if the gene is not cloned!
        }
    }
}

static utype measure_distance_stages(const struct Pipeline *pipeline, utype gene_id, utype stage_id_a, utype stage_id_b)
{
    assert(stage_id_a < pipeline->stage_count && stage_id_b < pipeline->stage_count && stage_id_a < stage_id_b && gene_id < pipeline->max_gene);
    const struct Gene *current_gene = &pipeline->genes[gene_id];
    utype distance = 0;
    for_loop(i, pipeline->atom_count)
    {
        utype affiliated_group_a = pipeline->stage_infos[stage_id_a].atom_affiliations[i];
        utype affiliated_group_b = pipeline->stage_infos[stage_id_b].atom_affiliations[i];
        const struct vec2u *start_relative_pos = &current_gene->atom_relative_pos[stage_id_a][i];
        const struct vec2u *target_relative_pos = &current_gene->atom_relative_pos[stage_id_b][i];
        const struct vec2u *start_group_pos = &current_gene->group_pos[stage_id_a][affiliated_group_a];
        const struct vec2u *target_group_pos = &current_gene->group_pos[stage_id_b][affiliated_group_b];
        utype start_pos_x = start_relative_pos->x + start_group_pos->x * pipeline->grid_config.ria_size.x;
        utype start_pos_y = start_relative_pos->y + start_group_pos->y * pipeline->grid_config.ria_size.y;
        utype target_pos_x = target_relative_pos->x + target_group_pos->x * pipeline->grid_config.ria_size.x;
        utype target_pos_y = target_relative_pos->y + target_group_pos->y * pipeline->grid_config.ria_size.y;
        long long gradx = (long long)target_pos_x - (long long)start_pos_x;
        long long grady = (long long)target_pos_y - (long long)start_pos_y;
        distance += (utype)(llabs(gradx) + llabs(grady));
    }
    return distance;
}

void measure_all_scores(struct Pipeline *pipeline)
{
    for_loop(i, pipeline->max_gene)
    {
        //printf("test %lu\n", i);
        if (pipeline->genes[i].is_alive == false || pipeline->genes[i].score != INVALID_UTYPE) continue; // skip dead and "already measured" genes
        pipeline->genes[i].score = 0;
        //printf("score %u\n", pipeline->genes[i].score);
        //show_gene_to_console(pipeline, i, 0);
        for_loop(j, pipeline->stage_count - 1)
        {
            pipeline->genes[i].score += measure_distance_stages(pipeline, i, j, j + 1);
            //printf("inc %u\n", pipeline->genes[i].score);
        }
        //printf("%lu = %u\n", i, pipeline->genes[i].score);
    }
}

static utype get_survivor_with_max_score(const struct Pipeline *pipeline)
{
    utype max_survivor_id = 0;
    for(utype i = 1; i < pipeline->gene_alive_count; ++i)
    {
        if (pipeline->genes[pipeline->survivor_ids[i]].score > pipeline->genes[pipeline->survivor_ids[max_survivor_id]].score)
        {
            max_survivor_id = i;
        }
    }
    return max_survivor_id;
}

utype darwin(struct Pipeline *pipeline, bool show_generation_stats)
{
    /// debug
    //utype min_pre = INVALID_UTYPE, max_pre = 0;
    //for_loop(i, pipeline->max_gene)
    //{
    //    if (pipeline->genes[i].score < min_pre) min_pre = pipeline->genes[i].score;
    //    if (pipeline->genes[i].score > max_pre) max_pre = pipeline->genes[i].score;
    //}
    //printf("pre min = %u | pre max = %u\n", min_pre, max_pre);
    /// fill survivor list with first genes
    utype survivor_i = 0, gene_i = 0;
    while (survivor_i < pipeline->gene_alive_count && gene_i < pipeline->max_gene)
    {
        pipeline->genes[gene_i].is_survivor = false; // reseting survivor flag
        if (pipeline->genes[gene_i].is_alive == false)
        {
            gene_i += 1;
            continue;
        }
        pipeline->survivor_ids[survivor_i] = gene_i;
        gene_i += 1;
        survivor_i += 1;
    }
    /// bug?
    if (gene_i >= pipeline->max_gene - 1)
    {
        printf("not enough alive genes... bug");
        exit(1);
    }
    /// pseudo-sorting
    utype max_survivor_id = get_survivor_with_max_score(pipeline);
    for (utype i = gene_i; i < pipeline->max_gene; ++i)
    {
        pipeline->genes[i].is_survivor = false; // reseting survivor flag
        if (pipeline->genes[i].is_alive == false) continue;
        if (pipeline->genes[i].score < pipeline->genes[pipeline->survivor_ids[max_survivor_id]].score)
        {
            pipeline->survivor_ids[max_survivor_id] = i;
            max_survivor_id = get_survivor_with_max_score(pipeline);
        }
    }
    /// turn flags for the survivor and killed genes
    utype min_id = INVALID_UTYPE, min_score = INVALID_UTYPE, max_score = 0, sum_score = 0;
    for_loop(i, pipeline->gene_alive_count)
    {
        struct Gene *current_gene = &pipeline->genes[pipeline->survivor_ids[i]];
        current_gene->is_survivor = true;
        current_gene->is_alive = true;
        if (min_score > current_gene->score)
        {
            min_score = current_gene->score;
            min_id = pipeline->survivor_ids[i];
        }
        if (max_score < current_gene->score) max_score = current_gene->score;
        sum_score += current_gene->score;
    }
    if (show_generation_stats == true)
    {
        float average_score = (float)sum_score / (float)pipeline->gene_alive_count;
        printf("min = %u | max = %u | average = %.2f\n", min_score, max_score, average_score);
    }
    for_loop(i, pipeline->max_gene)
    {
        if (pipeline->genes[i].is_survivor == false)
        {
            pipeline->genes[i].is_alive = false;
        }
    }
    return min_id;
}

void show_gene_to_console(const struct Pipeline *pipeline, utype gene_id, utype stage_id)
{
    show_gene_state_to_console(pipeline, gene_id);
    assert(stage_id < pipeline->stage_count);
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

void show_gene_state_to_console(const struct Pipeline *pipeline, utype gene_id)
{
    assert(gene_id < pipeline->max_gene);
    printf("Gene %u Score %u Survivor %i Alive %i\n", gene_id, pipeline->genes[gene_id].score, pipeline->genes[gene_id].is_survivor, pipeline->genes[gene_id].is_alive);
}

void show_all_gene_states_to_console(const struct Pipeline *pipeline)
{
    printf("All genes...\n");
    for_loop(i, pipeline->max_gene)
    {
        show_gene_state_to_console(pipeline, i);
    }
}
