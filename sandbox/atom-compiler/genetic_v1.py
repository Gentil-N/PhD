from enum import Enum
from typing import List
from typing import Self
from typing import Dict
import copy
import numpy as np
import random
import math
from vec2 import Vec2
from itertools import chain
import matplotlib.pyplot as plt
import networkx as nx
from sklearn.manifold import MDS



RIA_SIDE_SIZE = 2
GRID_SLIGHT_SHUFFLE_DISTANCE = 2
RIA_SLIGHT_SHUFFLE_DISTANCE = 2



def rlf(graph: nx.Graph) -> int:
    partition = []
    current_set = []
    node_list = list(graph.nodes)
    #print("node list ", node_list)
    while len(node_list) != 0:
        current_set.append(node_list[0])
        for node in node_list:
            if len(graph.adj[node]) > len(graph.adj[current_set[0]]):
                current_set[0] = node
        node_list.remove(current_set[0])
        independent_node_tuples = []
        for node in node_list:
            if not(current_set[0] in graph.adj[node]):
                independent_node_tuples.append((node, 0, 0))
        while len(independent_node_tuples) != 0:
            for node_tuple in independent_node_tuples:
                for neighbor in graph.adj[node_tuple[0]]:
                    if neighbor in current_set:
                        node_tuple = (node_tuple[0], node_tuple[1] + 1, node_tuple[2]) # improve perf
                    else:
                        node_tuple = (node_tuple[0], node_tuple[1], node_tuple[2] + 1) # improve perf
            max_list = [independent_node_tuples[0]]
            for node_tuple in independent_node_tuples:
                if node_tuple[1] < max_list[0][1]:
                    continue
                elif node_tuple[1] > max_list[0][1]:
                    max_list.clear()
                max_list.append(node_tuple)
            selected_node_tuple = max_list[0]
            for node_tuple in max_list:
                if node_tuple[2] < selected_node_tuple[2]:
                    selected_node_tuple = node_tuple
            independent_node_tuples.remove(selected_node_tuple)
            current_set.append(selected_node_tuple[0])
            independent_node_tuples = [node_tuple for node_tuple in independent_node_tuples if not selected_node_tuple[0] in graph.adj[node_tuple[0]]]
        graph.remove_nodes_from(current_set)
        partition.append(current_set.copy())
        current_set.clear()
        node_list = list(graph.nodes)
    #print(partition)
    return len(partition)



#print("start")
#g = nx.Graph()
#g.add_node(0)
#g.add_node(1)
#g.add_node(2)
#g.add_node(3)
#g.add_node(4)
#g.add_node(5)
#g.add_edge(0, 3)
#g.add_edge(0, 4)
#g.add_edge(0, 5)
#g.add_edge(1, 3)
#g.add_edge(1, 4)
#g.add_edge(1, 5)
#g.add_edge(2, 3)
#g.add_edge(2, 4)
#g.add_edge(2, 5)
#print(rlf(g))
#exit()



def rand2Dexcept(current: Vec2, shift: int, upper_bound_x: int, lower_bound_x: int, upper_bound_y: int, lower_bound_y: int) -> Vec2:
    randx = random.choice(range(max(lower_bound_x, current.x - shift), min(upper_bound_x, current.x + shift) + 1))
    if randx == current.x:
        y_list = list(range(max(lower_bound_y, current.y - shift), min(upper_bound_y, current.y + shift) + 1))
        y_list.remove(current.y)
        randy = random.choice(y_list)
        return Vec2(randx, randy)
    else:
        y_range = range(max(lower_bound_y, current.y - shift), min(upper_bound_y, current.y + shift) + 1)
        return Vec2(randx, random.choice(y_range))



class Atom:

    def __init__(self, id: int = -1, pos: Vec2 = Vec2(-1, -1)) -> None:
        self.id: int = id
        self.pos: Vec2 = pos

    def is_valid(self):
        return self.id == -1



class Ria:

    def __init__(self, pos: Vec2 = Vec2(-1, -1)) -> None:
        self.pos = pos
        self.atoms: List[Atom] = []

    def is_pos_occupied(self, pos: Vec2) -> tuple[bool, Atom]:
        for atom in self.atoms:
            if atom.pos.x == pos.x and atom.pos.y == pos.y:
                return (True, atom)
        return (False, Atom())

    def has_atom(self, atom_id: int) -> tuple[bool, Atom]:
        for atom in self.atoms:
            if atom.id == atom_id:
                return (True, atom)
        return (False, Atom())

    def add_atom(self, x: int, y: int, id: int):
        pos = Vec2(x, y)
        assert self.is_pos_occupied(pos)[0] == False
        self.atoms.append(Atom(id, pos))

    def switch_atoms(self, x_a: int, y_a: int, x_b: int, y_b: int):
        res_a = self.is_pos_occupied(Vec2(x_a, y_a))
        assert res_a[0] == True
        res_b = self.is_pos_occupied(Vec2(x_b, y_b))
        if res_b[0] == True:
            res_a[1].pos.x, res_b[1].pos.x = res_b[1].pos.x, res_a[1].pos.x
            res_a[1].pos.y, res_b[1].pos.y = res_b[1].pos.y, res_a[1].pos.y
        else:
            res_a[1].pos.x = x_b
            res_a[1].pos.y = y_b

    def apply_slight_shuffle(self):
        atom_in_use = random.choice(self.atoms)
        new_pos = rand2Dexcept(atom_in_use.pos, RIA_SLIGHT_SHUFFLE_DISTANCE, RIA_SIDE_SIZE - 1, 0, RIA_SIDE_SIZE - 1, 0)
        #print(atom_in_use.pos.x, atom_in_use.pos.y, new_pos.x, new_pos.y)
        self.switch_atoms(atom_in_use.pos.x, atom_in_use.pos.y, new_pos.x, new_pos.y)



class Grid:

    def __init__(self, ria_count: int) -> None:
        self.rias: List[List] = [[None for i in range(ria_count)] for j in range(ria_count)]
        self.ria_collected: List[Ria] = []

    def len(self):
        return len(self.rias)

    def lenlen(self):
        return len(self.rias[0])

    def is_ria_occupied(self, x: int, y: int) -> tuple[bool, Ria]:
        assert x >= 0 and x < self.len() and y >= 0 and y < self.lenlen()
        return (self.rias[x][y] != None, self.rias[x][y])

    def put_ria(self, ria: Ria):
        assert ria.pos.x >= 0 and ria.pos.x < self.len() and ria.pos.y >= 0 and ria.pos.y < self.lenlen()
        assert self.is_ria_occupied(ria.pos.x, ria.pos.y)[0] == False
        self.rias[ria.pos.x][ria.pos.y] = ria

    def get_atom_from_absolute_pos(self, x: int, y: int):
        assert x >= 0 and x < self.len() * RIA_SIDE_SIZE and y >= 0 and y < self.lenlen() * RIA_SIDE_SIZE
        res = self.is_ria_occupied(math.floor(x / RIA_SIDE_SIZE), math.floor(y / RIA_SIDE_SIZE))
        if res[0] == False:
            return None
        subres = res[1].is_pos_occupied(Vec2(x - res[1].pos.x * RIA_SIDE_SIZE,  y - res[1].pos.y * RIA_SIDE_SIZE))
        if subres[0] == False:
            return None
        return subres[1]

    def find_atom(self, atom_id: int) -> tuple[Ria, Atom]:
        for i in range(self.len()):
            for j in range(self.lenlen()):
                res = self.is_ria_occupied(i, j)
                if res[0] == False:
                    continue
                subres = res[1].has_atom(atom_id)
                if subres[0] == False:
                    continue
                return (res[1], subres[1])
        raise Exception("unable to find atom " + str(atom_id))

    def collect_rias(self):
        self.ria_collected.clear()
        for i in range(self.len()):
            for j in range(self.lenlen()):
                if self.rias[i][j] != None:
                    self.ria_collected.append(self.rias[i][j])

    def switch_rias(self, x_a: int, y_a: int, x_b: int, y_b: int):
        res_a = self.is_ria_occupied(x_a, y_a)
        assert res_a[0] == True
        res_b = self.is_ria_occupied(x_b, y_b)
        if res_b[0] == True:
            res_a[1].pos.x, res_b[1].pos.x = res_b[1].pos.x, res_a[1].pos.x
            res_a[1].pos.y, res_b[1].pos.y = res_b[1].pos.y, res_a[1].pos.y
        else:
            res_a[1].pos.x = x_b
            res_a[1].pos.y = y_b
        self.rias[x_a][y_a], self.rias[x_b][y_b] = res_b[1], res_a[1]

    def apply_slight_shuffle_on_grid(self):
        ria_in_use: Ria = random.choice(self.ria_collected)
        new_pos = rand2Dexcept(ria_in_use.pos, GRID_SLIGHT_SHUFFLE_DISTANCE, self.len() - 1, 0, self.lenlen() - 1, 0)
        self.switch_rias(ria_in_use.pos.x, ria_in_use.pos.y, new_pos.x, new_pos.y)

    def apply_slight_shuffle_on_random_ria(self):
        ria_in_use: Ria = random.choice(self.ria_collected)
        ria_in_use.apply_slight_shuffle()

    def print(self):
        res = ""
        for i in range(self.len() * RIA_SIDE_SIZE):
            for j in range(self.lenlen() * RIA_SIDE_SIZE):
                atom = self.get_atom_from_absolute_pos(i, j)
                if atom == None:
                    res += "."
                else:
                    res += str(atom.id)
                res += "\t"
            res += "\n"
        print(res)

    def encode(self) -> List[int]:
        assert len(self.ria_collected) > 0
        encoded: List[int] = []
        for ria in self.ria_collected:
            for atom in ria.atoms:
                encoded.extend([ria.pos.x * RIA_SIDE_SIZE + atom.pos.x, ria.pos.y * RIA_SIDE_SIZE + atom.pos.y])
        return encoded



def estimate_transfer_cost(grid_a: Grid, grid_b: Grid) -> int:
    cost: int = 0
    for i in range(grid_a.len()):
        for j in range(grid_a.lenlen()):
            res = grid_a.is_ria_occupied(i, j)
            if res[0] == False:
                continue
            for atom_a in res[1].atoms:
                ria_atom_b = grid_b.find_atom(atom_a.id)
                temp_cost = abs((atom_a.pos.x + res[1].pos.x * RIA_SIDE_SIZE) - (ria_atom_b[1].pos.x + ria_atom_b[0].pos.x * RIA_SIDE_SIZE)) + abs((atom_a.pos.y + res[1].pos.y * RIA_SIDE_SIZE) - (ria_atom_b[1].pos.y + ria_atom_b[0].pos.y * RIA_SIDE_SIZE))
                cost += temp_cost
    return cost



def are_atom_moves_compatible(atom_a_pos_a: Vec2, atom_a_pos_b: Vec2, atom_b_pos_a: Vec2, atom_b_pos_b: Vec2):
    if atom_a_pos_a.x < atom_b_pos_a.x and atom_a_pos_b.x < atom_b_pos_b.x:
        if atom_a_pos_a.y < atom_b_pos_a.y and atom_a_pos_b.y < atom_b_pos_b.y:
            return True
        if atom_a_pos_a.y > atom_b_pos_a.y and atom_a_pos_b.y > atom_b_pos_b.y:
            return True
        if atom_a_pos_a.y == atom_b_pos_a.y and atom_a_pos_b.y == atom_b_pos_b.y:
            return True
    if atom_a_pos_a.x > atom_b_pos_a.x and atom_a_pos_b.x > atom_b_pos_b.x:
        if atom_a_pos_a.y < atom_b_pos_a.y and atom_a_pos_b.y < atom_b_pos_b.y:
            return True
        if atom_a_pos_a.y > atom_b_pos_a.y and atom_a_pos_b.y > atom_b_pos_b.y:
            return True
        if atom_a_pos_a.y == atom_b_pos_a.y and atom_a_pos_b.y == atom_b_pos_b.y:
            return True
    if atom_a_pos_a.x == atom_b_pos_a.x and atom_a_pos_b.x == atom_b_pos_b.x:
        if atom_a_pos_a.y < atom_b_pos_a.y and atom_a_pos_b.y < atom_b_pos_b.y:
            return True
        if atom_a_pos_a.y > atom_b_pos_a.y and atom_a_pos_b.y > atom_b_pos_b.y:
            return True
        if atom_a_pos_a.y == atom_b_pos_a.y and atom_a_pos_b.y == atom_b_pos_b.y:
            raise Exception("unexpected bug: both atoms at same place? ")
    return False



def fill_compatibility_group(id: int, atom_pos_a: Vec2, atom_pos_b: Vec2, groups: List[List[tuple[Vec2, Vec2]]], id_groups: List[List[int]]):
    for i in range(len(groups)):
        compatibility_flag = True
        for elem in groups[i]:
            if not are_atom_moves_compatible(atom_pos_a, atom_pos_b, elem[0], elem[1]):
                compatibility_flag = False
                break
        if compatibility_flag:
            groups[i].append((atom_pos_a, atom_pos_b))
            id_groups[i].append(id)
            return # the atom is placed in a group: the job is finished
    groups.append([(atom_pos_a, atom_pos_b)]) # if no suitable group found: append a new group with the specific atom move
    id_groups.append([id])



def estimate_group_count(grid_a: Grid, grid_b: Grid) -> tuple[int, List[List[int]]]:
    groups: List[List[tuple[Vec2, Vec2]]] = []
    id_groups: List[List[int]] = []
    for ria in grid_a.ria_collected:
        for atom_grid_a in ria.atoms:
            res_grid_b = grid_b.find_atom(atom_grid_a.id)
            fill_compatibility_group(atom_grid_a.id, Vec2(atom_grid_a.pos.x + ria.pos.x * RIA_SIDE_SIZE, atom_grid_a.pos.y + ria.pos.y * RIA_SIDE_SIZE), Vec2(res_grid_b[1].pos.x + res_grid_b[0].pos.x * RIA_SIDE_SIZE, res_grid_b[1].pos.y + res_grid_b[0].pos.y * RIA_SIDE_SIZE), groups, id_groups)
    return (len(groups), id_groups)



def estimate_gradient_disparity(grid_a: Grid, grid_b: Grid) -> float:
    assert grid_a.len() == grid_b.len() and grid_a.lenlen() == grid_b.lenlen()
    x_disparity: List[int] = [0] * ((grid_a.len() * RIA_SIDE_SIZE - 1) * 2 + 1)
    y_disparity: List[int] = [0] * ((grid_a.lenlen() * RIA_SIDE_SIZE - 1) * 2 + 1)
    for ria_a in grid_a.ria_collected:
        for atom_a in ria_a.atoms:
            atom_b = grid_b.find_atom(atom_a.id)
            origin_x = ria_a.pos.x * RIA_SIDE_SIZE + atom_a.pos.x
            origin_y = ria_a.pos.y * RIA_SIDE_SIZE + atom_a.pos.y
            dir_x = atom_b[0].pos.x * RIA_SIDE_SIZE + atom_b[1].pos.x - origin_x
            dir_y = atom_b[0].pos.y * RIA_SIDE_SIZE + atom_b[1].pos.y - origin_y
            #print(dir_x, dir_y)
            x_disparity[dir_x + grid_a.len() * RIA_SIDE_SIZE - 1] += 1
            y_disparity[dir_y + grid_a.lenlen() * RIA_SIDE_SIZE - 1] += 1
    total_x_disparity = 0
    for num in x_disparity:
        if num > 0:
            total_x_disparity += 1
    total_y_disparity = 0
    for num in y_disparity:
        if num > 0:
            total_y_disparity += 1
    return (float(total_x_disparity) + float(total_y_disparity)) / 2.0



def estimate_coloring_rlf(grid_a: Grid, grid_b: Grid) -> int:
    assert grid_a.len() == grid_b.len() and grid_a.lenlen() == grid_b.lenlen()
    graph = nx.Graph()
    moves: List[tuple[Ria, Atom, Ria, Atom]] = []
    # Build nodes
    for ria_a in grid_a.ria_collected:
        for atom_a in ria_a.atoms:
            atom_b = grid_b.find_atom(atom_a.id)
            graph.add_node(atom_a.id)
            moves.append((ria_a, atom_a, atom_b[0], atom_b[1]))
    # Build edges
    for i in range(len(moves)):
        for j in range(i + 1, len(moves)):
            atom_a_pos_a = Vec2(moves[i][1].pos.x + moves[i][0].pos.x * RIA_SIDE_SIZE, moves[i][1].pos.y + moves[i][0].pos.y * RIA_SIDE_SIZE)
            atom_a_pos_b = Vec2(moves[i][3].pos.x + moves[i][2].pos.x * RIA_SIDE_SIZE, moves[i][3].pos.y + moves[i][2].pos.y * RIA_SIDE_SIZE)
            atom_b_pos_a = Vec2(moves[j][1].pos.x + moves[j][0].pos.x * RIA_SIDE_SIZE, moves[j][1].pos.y + moves[j][0].pos.y * RIA_SIDE_SIZE)
            atom_b_pos_b = Vec2(moves[j][3].pos.x + moves[j][2].pos.x * RIA_SIDE_SIZE, moves[j][3].pos.y + moves[j][2].pos.y * RIA_SIDE_SIZE)
            if not are_atom_moves_compatible(atom_a_pos_a, atom_a_pos_b, atom_b_pos_a, atom_b_pos_b):
                graph.add_edge(moves[i][1].id, moves[j][1].id)
    #nx.draw(graph)
    #plt.show()
    return rlf(graph)



def generate_grid(ria_count: int, groups: List[List[int]]) -> Grid:
    assert len(groups) <= ria_count**2
    grid: Grid = Grid(ria_count)
    for i in range(len(groups)):
        assert len(groups[i]) <= RIA_SIDE_SIZE**2
        ria: Ria = Ria(Vec2(math.floor(i / ria_count), i % ria_count))
        for j in range(len(groups[i])):
            ria.add_atom(math.floor(j / RIA_SIDE_SIZE), j % RIA_SIDE_SIZE, groups[i][j])
        grid.put_ria(ria)
    grid.collect_rias()
    return grid



class Family:

    def __init__(self, ria_count: int, stages: List[List[List[int]]]) -> None:
        assert len(stages) >= 2
        self.family: List[Grid] = []
        self.encoded = False
        for stage in stages:
            self.family.append(generate_grid(ria_count, stage))

    def shuffle(self, ria_pass_count: int, atom_pass_count: int):
        for grid in self.family:
            for i in range(ria_pass_count):
                grid.apply_slight_shuffle_on_grid()
            for j in range(atom_pass_count):
                grid.apply_slight_shuffle_on_random_ria()

    def copy(self) -> Self:
        cp = copy.deepcopy(self)
        cp.encoded = False
        return cp

    def estimate_total_transfer_cost(self) -> int:
        total_cost = 0
        for i in range(len(self.family) - 1):
            total_cost += estimate_transfer_cost(self.family[i], self.family[i + 1])
        return total_cost

    def estimate_total_group_count(self) -> int:
        total_group_count = 0
        for i in range(len(self.family) - 1):
            total_group_count += estimate_group_count(self.family[i], self.family[i + 1])[0]
        return total_group_count

    def estimate_total_gradient_disparity(self) -> int:
        total_disparity = 0
        for i in range(len(self.family) - 1):
            total_disparity += estimate_gradient_disparity(self.family[i], self.family[i + 1])
        return round(total_disparity)

    def estimate_total_gradient_disparity_for_transfer(self, transfer_id: int) -> float:
        assert transfer_id < len(self.family) - 1
        return estimate_gradient_disparity(self.family[transfer_id], self.family[transfer_id + 1])

    def estimate_total_coloring_rlf(self) -> int:
        total_coloring = 0
        for i in range(len(self.family) - 1):
            total_coloring += estimate_coloring_rlf(self.family[i], self.family[i + 1])
        return round(total_coloring)

    def plot_gradient_field(self, transfer_id: int, type: int):
        assert transfer_id < len(self.family) - 1
        stage_a_id = transfer_id
        stage_b_id = transfer_id + 1
        for ria_a in self.family[stage_a_id].ria_collected:
            for atom_a in ria_a.atoms:
                atom_b = self.family[stage_b_id].find_atom(atom_a.id)
                origin_x = ria_a.pos.x * RIA_SIDE_SIZE + atom_a.pos.x
                origin_y = ria_a.pos.y * RIA_SIDE_SIZE + atom_a.pos.y
                dest_x = atom_b[0].pos.x * RIA_SIDE_SIZE + atom_b[1].pos.x
                dest_y = atom_b[0].pos.y * RIA_SIDE_SIZE + atom_b[1].pos.y
                dir_x = dest_x - origin_x
                dir_y = dest_y - origin_y
                if type == 0:
                    plt.arrow(0, 0, dir_x, dir_y, width=0.03)
                elif type == 1:
                    plt.arrow(origin_x, origin_y, dir_x, dir_y, width=0.03)
                else:
                    raise Exception("invalid plot type")
        if type == 0:
            plt.ylim(-self.family[stage_a_id].len() * RIA_SIDE_SIZE, self.family[stage_a_id].len() * RIA_SIDE_SIZE)
            plt.xlim(-self.family[stage_a_id].lenlen() * RIA_SIDE_SIZE, self.family[stage_a_id].lenlen() * RIA_SIDE_SIZE)
        elif type == 1:
            plt.ylim(-1, self.family[stage_a_id].len() * RIA_SIDE_SIZE + 1)
            plt.xlim(-1, self.family[stage_a_id].lenlen() * RIA_SIDE_SIZE + 1)
        plt.show()

    def print(self):
        for i in range(len(self.family)):
            print("Stage ", i)
            self.family[i].print()

    def encode_two_stages(self, first_stage_index: int = 0) -> List[int]:
        assert first_stage_index < len(self.family) - 1
        encoded = self.family[first_stage_index].encode()
        encoded.extend(self.family[first_stage_index + 1].encode())
        return encoded

    def encode_all_stages(self) -> List[int]:
        assert self.encoded == False
        encoded: List[int] = []
        for grid in self.family:
            encoded.extend(grid.encode())
        self.encoded = True
        return encoded

    def is_encoded(self) -> bool:
        return self.encoded



def genesis(family_count: int, ria_count: int, ria_pass_count: int, atom_pass_count: int, stages: List[List[List[int]]]) -> List[Family]:
    families: List[Family] = []
    for i in range(family_count):
        families.append(Family(ria_count, stages))
        families[-1].shuffle(ria_pass_count, atom_pass_count)
    return families



def grow_population(population: List[Family], children_per_family: int, ria_pass_count: int, atom_pass_count: int):
    original_length = len(population)
    for i in range(original_length):
        for j in range(children_per_family):
            new_family = population[i].copy()
            new_family.shuffle(ria_pass_count, atom_pass_count)
            population.append(new_family)



def get_sort_key(family_cost: tuple[Family, int]) -> int:
    return family_cost[1]



def kill_worst_families(generation_id: int, population: List[Family], survivor_count: int, selection_rule: int, encoded, z_measurement, generation_id_colors) -> tuple[int, int, float, List[Family]]:
    assert len(population) >= 1 and survivor_count <= len(population)
    pop_cost: List[tuple[Family, int]] = []
    for family in population:
        if selection_rule == 0:
            pop_cost.append((family, family.estimate_total_transfer_cost()))
        elif selection_rule == 1:
            pop_cost.append((family, family.estimate_total_group_count()))
        elif selection_rule == 2:
            pop_cost.append((family, family.estimate_total_gradient_disparity()))
        elif selection_rule == 3:
            pop_cost.append((family, family.estimate_total_coloring_rlf()))
        else:
            raise Exception("selection rule not valid")
    pop_cost.sort(key=get_sort_key)
    survivors: List[Family] = []
    average: float = 0.0
    for i in range(min(survivor_count, len(pop_cost))):
        survivors.append(pop_cost[i][0])
        average += float(pop_cost[i][1])
    average /= len(survivors)
    if encoded != None and z_measurements != None and generation_id_colors != None:
        #print("encoding")
        for i in range(len(survivors)):
            if not survivors[i].is_encoded():
                #print("new one")
                encoded.append(survivors[i].encode_all_stages())
                z_measurement.append(pop_cost[i][1])
                generation_id_colors.append(generation_id)
    return (pop_cost[0][1], pop_cost[len(survivors) - 1][1], average, survivors)



def process_generation(generation_id: int, population: List[Family], children_per_family: int, ria_pass_count: int, atom_pass_count: int, selection_rule: int, encoded, z_measurement, generation_id_colors) -> tuple[int, int, float, List[Family]]:
    survivor_count: int = len(population)
    grow_population(population, children_per_family, ria_pass_count, atom_pass_count)
    return kill_worst_families(generation_id, population, survivor_count, selection_rule, encoded, z_measurement, generation_id_colors)



def darwin(genesis_population: List[Family], children_per_family: int, ria_pass_count: int, atom_pass_count: int, generation_limit : int, generation_shift: int, selection_rule: int, encoded, z_measurement, generation_id_colors) -> List[Family]:
    population = copy.deepcopy(genesis_population)
    for i in range(generation_limit):
        res = process_generation(generation_shift + i, population, children_per_family, ria_pass_count, atom_pass_count, selection_rule, encoded, z_measurement, generation_id_colors)
        print("Generation ", i, ": min = ", res[0], " | max = ", res[1], " | average = ", res[2])
        population = res[3]
    print("Best family found is...")
    #print("Estimated cost = ", population[0].estimate_total_transfer_cost())
    total_goup_count = 0
    for i in range(len(population[0].family) - 1):
        res = estimate_group_count(population[0].family[i], population[0].family[i + 1])
        total_goup_count += res[0]
        print("Group ", i, res[1])
    print("Estimated total group count = ", total_goup_count)
    population[0].print()
    return population



#def darwin_randpass(genesis_population: List[Family], children_per_family: int, ria_pass_count_limit: int, atom_pass_count_limit: int, generation_limit : int) -> List[Family]:
#    population = copy.deepcopy(genesis_population)
#    for i in range(generation_limit):
#        res = process_generation(population, children_per_family, random.choice(range(ria_pass_count_limit + 1)), random.choice(range(atom_pass_count_limit + 1)), 0)
#        print("Generation ", i, ": min = ", res[0], " | max = ", res[1], " | average = ", res[2])
#        population = res[3]
#    print("Best family found is...")
#    #print("Estimated cost = ", population[0].estimate_total_transfer_cost())
#    total_goup_count = 0
#    for i in range(len(population[0].family) - 1):
#        res = estimate_group_count(population[0].family[i], population[0].family[i + 1])
#        total_goup_count += res[0]
#        print("Group ", i, res[1])
#    print("Estimated total group count = ", total_goup_count)
#    population[0].print()
#    return population



def plot_mds(vectors: List[List[int]], z_measurements: List[int], colors: List[int]):
    print("processing mds")
    mds = MDS(n_components=2, metric=True)
    result = mds.fit_transform(vectors)
    print(mds.stress_)
    xy_reduced: np.ndarray = np.asarray(result)
    #print(xy_reduced.shape)
    print("plotting...")
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    #print(z_measurements)
    #print(xy_reduced)
    ax.scatter(xy_reduced[:, 0], xy_reduced[:, 1], z_measurements, c=colors, cmap="rainbow")
    plt.show()



# TEST



#ria_a_0 = Ria(Vec2(0, 0))
#ria_a_0.add_atom(0, 0, 0)
#ria_a_0.add_atom(0, 1, 1)
#ria_a_0.add_atom(1, 0, 2)
#ria_a_0.add_atom(1, 1, 3)
#grid_a = Grid(4)
#grid_a.put_ria(ria_a_0)
#grid_a.print()
#
#ria_b_0 = Ria(Vec2(0, 0))
#ria_b_0.add_atom(0, 0, 0)
#ria_b_0.add_atom(0, 1, 1)
#ria_b_0.add_atom(1, 0, 2)
#ria_b_0.add_atom(1, 1, 3)
#grid_b = Grid(4)
#grid_b.put_ria(ria_b_0)
#grid_b.print()
#
#print(estimate_transfer_cost(grid_a, grid_b))

#RIA_COUNT = 4
#stage_0: List[List[int]] = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9], [10, 11]]
#stage_1: List[List[int]] = [[0, 10, 11], [1, 2, 5], [3, 4, 6], [7, 8, 9]]
#stage_2: List[List[int]] = [[1, 5, 8, 10], [0, 2, 3], [4, 6, 9], [7, 11]]
#genesis_population = genesis(50, RIA_COUNT, 20, 20, [stage_0, stage_1, stage_2])
#population = darwin(genesis_population, 2, 1, 0, 500)
#population = darwin(population, 2, 0, 1, 500)

#RIA_COUNT = 4
#stage_0: List[List[int]] = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15], [16, 17], [18, 19], [20, 21], [22, 23], [24], [25]]
#stage_1: List[List[int]] = [[0, 1, 4, 5], [2, 3, 6, 7], [8, 9, 12, 13], [10, 11, 14, 15], [16, 20, 24], [17, 18, 19], [21, 22], [23], [25]]
#stage_2: List[List[int]] = [[0, 25, 8, 9], [1, 4, 6, 13], [2, 3, 10, 11], [5, 7, 12, 14], [15, 16], [20, 22], [17, 23], [24], [18], [19], [21]]
#genesis_population = genesis(50, RIA_COUNT, 20, 20, [stage_0, stage_1, stage_2])
##population = darwin(genesis_population, 2, 1, 0, 5000)
##population = darwin(population, 2, 0, 1, 5000)
#population = darwin(genesis_population, 2, 1, 0, 500)
#population = darwin(population, 2, 0, 1, 500)

RIA_COUNT = 4
stage_0: List[List[int]] = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15], [16, 17], [18, 19], [20, 21], [22, 23], [24], [25], [26, 27, 28], [29, 30, 31], [32, 33, 34], [35]]
stage_1: List[List[int]] = [[0, 1, 4, 5], [2, 3, 6, 7], [8, 9, 12, 13], [10, 11, 14, 15], [16, 20, 24], [17, 18, 19], [21, 22], [23], [25], [26, 29, 32, 35], [27, 30, 33], [28, 31, 34]]
stage_2: List[List[int]] = [[0, 25, 8, 9], [1, 4, 6, 13], [2, 3, 10, 11], [5, 7, 12, 14], [15, 16], [20, 22], [17, 23], [24], [18], [19], [21], [27, 32, 33, 35], [26, 30, 34], [28, 29, 31]]

encoded: List[List[int]] = []
z_measurements: List[int] = []
generation_id_colors: List[int] = []

def run_lap() -> int:
    population = genesis(2, RIA_COUNT, 20, 20, [stage_0, stage_1, stage_2])
    for family in population:
        encoded.append(family.encode_all_stages())
        z_measurements.append(family.estimate_total_transfer_cost())
        generation_id_colors.append(0)
    population = darwin(population, 1, 1, 0, 30, 0, 0, encoded, z_measurements, generation_id_colors)
    #population = darwin(population, 2, 0, 1, 200, 400, 0, encoded, z_measurements, generation_id_colors)
    return len(generation_id_colors)

colors: List[int] = []
last_size = 0
for i in range(1):
    current_size = run_lap()
    colors.extend([i] * (current_size - last_size))
    last_size = current_size

plot_mds(encoded, z_measurements, colors)
