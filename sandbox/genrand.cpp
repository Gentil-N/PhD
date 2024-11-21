#include "gencomp.h"

#include <cassert>
#include <vector>

#include <iostream>

#include "crand.h"

struct Group
{
    utype id;
    utype total_atom_count;
    utype current_atom_count;
};

utype generate_random_stage(utype atom_count, utype max_atom_per_group, utype *atom_affiliations)
{
    assert(atom_count > max_atom_per_group && max_atom_per_group > 0);
    crand_init();
    std::vector<Group> groups;
    utype total_atom_count = 0;
    utype id = 0;
    while (total_atom_count < atom_count)
    {
        utype current_atom_count_for_group = crand_generate(1, atom_count - total_atom_count >= max_atom_per_group ? max_atom_per_group : atom_count - total_atom_count);
        groups.push_back({id, current_atom_count_for_group, 0});
        //std::cout << current_atom_count_for_group << " ";
        total_atom_count += current_atom_count_for_group;
        id += 1;
    }
    //std::cout << std::endl;
    const utype group_count = (utype)groups.size();
    for(utype i = 0; i < atom_count; ++i)
    {
        utype group_index /* not group id! */ = crand_generate(0, (utype)groups.size() - 1);
        atom_affiliations[i] = groups[group_index].id;
        groups[group_index].current_atom_count += 1;
        if (groups[group_index].current_atom_count >= groups[group_index].total_atom_count) // if group full
        {
            groups.erase(groups.begin() + (size_t)group_index); // remove it from the list
        }
    }
    /// check sum
    //std::vector<utype> group_sum(group_count, 0);
    //for(size_t i = 0; i < atom_count; ++i)
    //{
    //    group_sum[atom_affiliations[i]] += 1;
    //    if (group_sum[atom_affiliations[i]] > max_atom_per_group)
    //    {
    //        std::cerr << "failed to generate random stage" << std::endl;
    //        exit(1);
    //    }
    //}
    //for (utype sum : group_sum)
    //{
    //    std::cout << sum << " ";
    //}
    //std::cout << std::endl;
    return group_count;
}
