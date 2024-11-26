#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "cutils.h"
#include "gencomp.h"

void trace_fn(void *address, const char *file, int line)
{
    printf("Memory leak (%p) in file %s at %d\n", address, filename(file), line);
}

/* reverse:  reverse string s in place */
void reverse(char s[])
{
    int i, j;
    char c;

    for (i = 0, j = strlen(s)-1; i<j; i++, j--) {
        c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}

/* itoa:  convert n to characters in s */
void itoa(int n, char s[])
{
    int i, sign;

    if ((sign = n) < 0)  /* record sign */
        n = -n;          /* make n positive */
    i = 0;
    do {       /* generate digits in reverse order */
        s[i++] = n % 10 + '0';   /* get next digit */
    } while ((n /= 10) > 0);     /* delete it */
    if (sign < 0)
        s[i++] = '-';
    s[i] = '\0';
    reverse(s);
}

void run_random(utype atom_count, int id)
{
    printf("Preparing random data\n");
    utype atom_affiliation_0[atom_count];
    utype group_count_0 = generate_random_stage(atom_count, 4, atom_affiliation_0);
    utype atom_affiliation_1[atom_count];
    utype group_count_1 = generate_random_stage(atom_count, 4, atom_affiliation_1);
    utype atom_affiliation_2[atom_count];
    utype group_count_2 = generate_random_stage(atom_count, 4, atom_affiliation_2);

    struct StageInfo stage_infos[] = {{group_count_0, NULL, atom_affiliation_0, NULL}, {group_count_1, NULL, atom_affiliation_1, NULL}, {group_count_2, NULL, atom_affiliation_2, NULL}};

    printf("Creating pipeline\n");
    struct GridConfig grid_config = {0};
    grid_config.ria_count = (struct vec2u){(utype)ceil(sqrtf((float)atom_count)), (utype)ceil(sqrtf((float)atom_count))};
    grid_config.ria_size = (struct vec2u){2, 2};
    struct Pipeline pipeline = {0};
    create_pipeline(&pipeline, &grid_config, atom_count, 1000, 2, 3, stage_infos);

    char *result = malloc(100 * sizeof(char));
    result[0] = '\0';
    char *str_id = malloc(100 * sizeof(char));
    itoa(id, str_id);
    //printf("%s\n", str_id);
    const char *base_log_name = "./log";
    strcat(result, base_log_name);
    strcat(result, str_id);
    struct LogFile *log_file = create_log_file(result);
    log_file_write_header(log_file, &pipeline);
    free(str_id);
    free(result);

    mutate_all_genes(&pipeline, 20, 20, 2, 2, false); // shuffle initial genes
    utype min_id = INVALID_UTYPE;
    for_loop(i, 500)
    {
        mutate_all_genes(&pipeline, 1, 0, 2, 2, true); // mutate
        measure_all_scores(&pipeline); // measure
        printf("Generation %lu ", i);
        min_id = darwin(&pipeline, true); // kill
    }
    for_loop(i, 500)
    {
        mutate_all_genes(&pipeline, 0, 1, 2, 2, true); // mutate
        measure_all_scores(&pipeline); // measure
        printf("Generation %lu ", i);
        min_id = darwin(&pipeline, true); // kill
    }

    //for_loop(j, pipeline.stage_count)
    //{
    //    show_gene_to_console(&pipeline, min_id, j);
    //}

    log_file_write_gene(log_file, &pipeline, min_id, "distance");

    destroy_log_file(log_file);
    destroy_pipeline(&pipeline);
}

int main()
{
    /*utype stage_0_group_count = 14;
    utype stage_0_atom_count_per_group[] = {4, 4, 4, 4, 2, 2, 2, 2, 1, 1, 3, 3, 3, 1};
    utype stage_0_atom_affiliation[] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12, 13};
    utype stage_0_data[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35};
    utype stage_1_group_count = 12;
    utype stage_1_atom_count_per_group[] = {4, 4, 4, 4, 3, 3, 2, 1, 1, 4, 3, 3};
    utype stage_1_atom_affiliation[] = {0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 3, 3, 2, 2, 3, 3, 4, 5, 5, 5, 4, 6, 6, 7, 4, 8, 9, 10, 11, 9, 10, 11, 9, 10, 11, 9};
    utype stage_1_data[] = {0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15, 16, 20, 24, 17, 18, 19, 21, 22, 23, 25, 26, 29, 32, 35, 27, 30, 33, 28, 31, 34};
    utype stage_2_group_count = 14;
    utype stage_2_atom_count_per_group[] = {4, 4, 4, 4, 2, 2, 2, 1, 1, 1, 1, 4, 3, 3};
    utype stage_2_atom_affiliation[] = {0, 1, 2, 2, 1, 3, 1, 3, 0, 0, 2, 2, 3, 1, 3, 4, 4, 6, 8, 9, 5, 10, 5, 6, 7, 0, 12, 11, 13, 13, 12, 13, 11, 11, 12, 11};
    utype stage_2_data[] = {0, 8, 9, 25, 1, 4, 6, 13, 2, 3, 10, 11, 5, 7, 12, 14, 15, 16, 20, 22, 17, 23, 24, 18, 19, 21, 27, 32, 33, 35, 26, 30, 34, 28, 29, 31};

    struct StageInfo stage_infos[] = {{stage_0_group_count, stage_0_atom_count_per_group, stage_0_atom_affiliation, stage_0_data}, {stage_1_group_count, stage_1_atom_count_per_group, stage_1_atom_affiliation, stage_1_data}, {stage_2_group_count, stage_2_atom_count_per_group, stage_2_atom_affiliation, stage_2_data}};

    struct GridConfig grid_config = {0};
    grid_config.ria_count = (struct vec2u){4, 4};
    grid_config.ria_size = (struct vec2u){2, 2};
    struct Pipeline pipeline = {0};
    create_pipeline(&pipeline, &grid_config, 36, 1000, 3, 3, stage_infos);

    struct LogFile *log_file = create_log_file("./log-test");
    log_file_write_header(log_file, &pipeline);

    mutate_all_genes(&pipeline, 20, 20, 2, 2, false); // shuffle initial genes
    utype min_id = INVALID_UTYPE;
    for_loop(i, 2000)
    {
        mutate_all_genes(&pipeline, 1, 0, 2, 2, true); // mutate
        measure_all_scores(&pipeline); // measure
        printf("Generation %lu ", i);
        min_id = darwin(&pipeline, true); // kill
    }
    for_loop(i, 2000)
    {
        mutate_all_genes(&pipeline, 0, 1, 2, 2, true); // mutate
        measure_all_scores(&pipeline); // measure
        printf("Generation %lu ", i);
        min_id = darwin(&pipeline, true); // kill
    }

    log_file_write_gene(log_file, &pipeline, min_id, "distance");

    //for_loop(j, pipeline.stage_count)
    //{
    //    show_gene_to_console(&pipeline, min_id, j);
    //}

    destroy_log_file(log_file);
    destroy_pipeline(&pipeline);*/



    /*utype atom_affiliation[36];
    utype group_count = generate_random_stage(36, 4, atom_affiliation);
    printf("group count = %u\n", group_count);
    for_loop(i, 36)
    {
        printf("%u ", atom_affiliation[i]);
    }
    printf("\n");
    exit(1);

    utype stage_0_group_count = 4;
    utype stage_0_atom_count_per_group[] = {4, 4, 4, 4};
    utype stage_0_atom_affiliation[] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
    utype stage_0_data[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    utype stage_1_group_count = 4;
    utype stage_1_atom_count_per_group[] = {4, 4, 4, 4};
    utype stage_1_atom_affiliation[] = {0, 1, 0, 1, 0, 1, 0, 1, 2, 3, 2, 3, 2, 3, 2, 3};
    utype stage_1_data[] = {0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    struct StageInfo stage_infos[] = {{stage_0_group_count, stage_0_atom_count_per_group, stage_0_atom_affiliation, stage_0_data}, {stage_1_group_count, stage_1_atom_count_per_group, stage_1_atom_affiliation, stage_1_data}};

    struct GridConfig grid_config = {0};
    grid_config.ria_count = (struct vec2u){6, 6};
    grid_config.ria_size = (struct vec2u){2, 2};
    struct Pipeline pipeline = {0};
    create_pipeline(&pipeline, &grid_config, 16, 20, 2, 2, stage_infos);

    //struct LogFile *log_file = create_log_file("./log-test");
    //log_file_write_header(log_file, &pipeline);

    mutate_all_genes(&pipeline, 20, 20, 2, 2, false); // shuffle initial genes
    utype min_id = INVALID_UTYPE;
    for_loop(i, 150)
    {
        mutate_all_genes(&pipeline, 1, 0, 2, 2, true); // mutate
        measure_all_scores(&pipeline); // measure
        printf("Generation %lu ", i);
        min_id = darwin(&pipeline, true); // kill
    }
    for_loop(i, 150)
    {
        mutate_all_genes(&pipeline, 0, 1, 2, 2, true); // mutate
        measure_all_scores(&pipeline); // measure
        printf("Generation %lu ", i);
        min_id = darwin(&pipeline, true); // kill
    }

    //log_file_write_gene(log_file, &pipeline, min_id, "distance");

    for_loop(j, pipeline.stage_count)
    {
        show_gene_to_console(&pipeline, min_id, j);
    }

    //destroy_log_file(log_file);
    destroy_pipeline(&pipeline);*/

    for(size_t i = 0; i < 1000; ++i)
    {
        run_random(9, (int)i);
    }
    tracker_trace(trace_fn);
}
