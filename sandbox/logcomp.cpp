#include "gencomp.h"

#include <cassert>
#include <string>
#include <fstream>

struct LogFile
{
    std::ofstream stream;
};

struct LogFile *create_log_file(const char *name)
{
    LogFile *log_file = new LogFile();
    log_file->stream.open(name);
    return log_file;
}

void destroy_log_file(struct LogFile *log_file)
{
    log_file->stream.close();
    delete log_file;
}

void log_file_write_header(struct LogFile *log_file, const struct Pipeline *pipeline)
{
    log_file->stream << "HEADER START" << std::endl;
    log_file->stream << "AtomComp Data File v" << ATOMCOMP_VERSION_MAJOR << "." << ATOMCOMP_VERSION_MINOR << std::endl;
    log_file->stream << "gridconfig=" << pipeline->grid_config.ria_count.x << "," << pipeline->grid_config.ria_count.y << "," << pipeline->grid_config.ria_size.x << "," << pipeline->grid_config.ria_size.y << std::endl;
    log_file->stream << "atomcount=" << pipeline->atom_count << std::endl;
    log_file->stream << "stagecount=" << pipeline->stage_count << std::endl;
    for(size_t i = 0; i < pipeline->stage_count; ++i)
    {
        //log_file->stream << "stageid=" << i << std::endl;
        log_file->stream << "groupcount" << i << "=" << pipeline->stage_infos[i].group_count << std::endl;
        log_file->stream << "atomaffiliations" << i << "=";
        for(size_t j = 0; j < pipeline->atom_count; ++j)
        {
            log_file->stream << pipeline->stage_infos[i].atom_affiliations[j] << (j == pipeline->atom_count - 1 ? "\n" : ",");
        }
    }
    log_file->stream << "HEADER END" << std::endl;
}

void log_file_write_gene(struct LogFile *log_file, const struct Pipeline *pipeline, utype gene_id, const char *selection_rule_name)
{
    assert(gene_id < pipeline->max_gene);
    Gene *current_gene = &pipeline->genes[gene_id];
    log_file->stream << "GENE START" << std::endl;
    log_file->stream << "geneid=" << gene_id << std::endl;
    log_file->stream << "selectionrule=" << selection_rule_name << std::endl;
    log_file->stream << "score=" << current_gene->score << std::endl;
    for(size_t i = 0; i < pipeline->stage_count; ++i)
    {
        //log_file->stream << "stageid=" << i << std::endl;
        for(size_t j = 0; j < pipeline->atom_count; ++j)
        {
            utype affiliated_group = pipeline->stage_infos[i].atom_affiliations[j];
            const struct vec2u *relative_pos = &current_gene->atom_relative_pos[i][j];
            const struct vec2u *group_pos = &current_gene->group_pos[i][affiliated_group];
            utype pos_x = relative_pos->x + group_pos->x * pipeline->grid_config.ria_size.x;
            utype pos_y = relative_pos->y + group_pos->y * pipeline->grid_config.ria_size.y;
            log_file->stream << i << "_" << j << "=" << pos_x << "," << pos_y << std::endl;
        }
    }
    log_file->stream << "GENE END" << std::endl;
}
