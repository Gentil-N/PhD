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

    def __init__(self, id: int = -1, pos: Vec2 = Vec2(-1, -1), ria_parent = None) -> None:
        self.id: int = id
        self.pos: Vec2 = pos
        self.gradient: Vec2 = Vec2(0.0, 0.0)
        self.ria_parent = None

    def is_valid(self):
        return self.id == -1



class Ria:

    def __init__(self, pos: Vec2 = Vec2(-1, -1)) -> None:
        self.pos = pos
        self.gradient: Vec2 = Vec2(0.0, 0.0)
        self.atoms: List[Atom] = []

    def is_pos_occupied(self, pos: Vec2) -> tuple[bool, Atom]:
        if pos.x < 0 or pos.y < 0 or pos.x >= RIA_SIDE_SIZE or pos.y >= RIA_SIDE_SIZE:
            return (False, Atom())
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
        self.atoms.append(Atom(id, pos, self))

    def switch_atoms(self, x_a: int, y_a: int, x_b: int, y_b: int):
        res_a = self.is_pos_occupied(Vec2(x_a, y_a))
        assert res_a[0] == True
        res_b = self.is_pos_occupied(Vec2(x_b, y_b))
        target_a = Vec2(res_a[1].pos.x + res_a[1].gradient.x, res_a[1].pos.y + res_a[1].gradient.y)
        if res_b[0] == True:
            target_b = Vec2(res_b[1].pos.x + res_b[1].gradient.x, res_b[1].pos.y + res_b[1].gradient.y)
            res_a[1].pos.x, res_b[1].pos.x = res_b[1].pos.x, res_a[1].pos.x
            res_a[1].pos.y, res_b[1].pos.y = res_b[1].pos.y, res_a[1].pos.y
            res_b[1].gradient.set_xy(target_b.x - res_b[1].pos.x, target_b.y - res_b[1].pos.y)
        else:
            res_a[1].pos.x = x_b
            res_a[1].pos.y = y_b
        res_a[1].gradient.set_xy(target_a.x - res_a[1].pos.x, target_a.y - res_a[1].pos.y)

    def apply_slight_shuffle(self):
        atom_in_use = random.choice(self.atoms)
        new_pos = rand2Dexcept(atom_in_use.pos, RIA_SLIGHT_SHUFFLE_DISTANCE, RIA_SIDE_SIZE - 1, 0, RIA_SIDE_SIZE - 1, 0)
        self.switch_atoms(atom_in_use.pos.x, atom_in_use.pos.y, new_pos.x, new_pos.y)

    def get_horizontal_best_swap(self, atom: Atom):
        left_coord = Vec2(atom.pos.x - 1, atom.pos.y)
        right_coord = Vec2(atom.pos.x + 1, atom.pos.y)
        left_tuple = self.is_pos_occupied(left_coord)
        right_tuple = self.is_pos_occupied(right_coord)
        to_choose: List[tuple[float, Vec2]] = []
        if left_tuple[0] == True: # left ria present
            if left_tuple[1] != None and left_tuple[1].gradient.x > atom.gradient.x: # if left ria filled and not ordered
                to_choose.append((abs(left_tuple[1].gradient.x - atom.gradient.x), left_coord))
            elif atom.gradient.x < 0 and abs(atom.gradient.x) > RIA_SIDE_SIZE / 2: # if not filled but present ria must move left
                to_choose.append((abs(atom.gradient.x), left_coord))
        if right_tuple[0] == True: # right ria present
            if right_tuple[1] != None and right_tuple[1].gradient.x < atom.gradient.x: # if right ria filled and not ordered
                to_choose.append((abs(right_tuple[1].gradient.x - atom.gradient.x), right_coord))
            elif atom.gradient.x > 0 and abs(atom.gradient.x) > RIA_SIDE_SIZE / 2: # if not filled but present ria must move right
                to_choose.append((abs(atom.gradient.x), right_coord))
        if len(to_choose) == 0: # no swap required
            return None
        if len(to_choose) == 1: # either the top or bottow ria is invalid or only top/bottom requires to be swaped
            return to_choose[0]
        if to_choose[0][0] > to_choose[1][0]:
            return to_choose[0]
        else:
            return to_choose[1]

    def get_vertical_best_swap(self, atom: Atom):
        bottom_coord = Vec2(atom.pos.x, atom.pos.y - 1)
        top_coord = Vec2(atom.pos.x, atom.pos.y + 1)
        bottom_tuple = self.is_pos_occupied(bottom_coord)
        top_tuple = self.is_pos_occupied(top_coord)
        to_choose: List[tuple[float, Vec2]] = []
        if bottom_tuple[0] == True: # bottom ria present
            if bottom_tuple[1] != None and bottom_tuple[1].gradient.y > atom.gradient.y: # if bottom ria filled and not ordered and diff grad greater than
                to_choose.append((abs(bottom_tuple[1].gradient.y - atom.gradient.y), bottom_coord))
            elif atom.gradient.y < 0 and abs(atom.gradient.y) > RIA_SIDE_SIZE / 2: # if not filled but the present ria must move bottom and gradient greater than atom count (side)
                to_choose.append((abs(atom.gradient.y), bottom_coord))
        if top_tuple[0] == True: # top ria present
            if top_tuple[1] != None and top_tuple[1].gradient.y < atom.gradient.y: # if top ria filled and not ordered
                to_choose.append((abs(top_tuple[1].gradient.y - atom.gradient.y), top_coord))
            elif atom.gradient.y > 0 and abs(atom.gradient.y) > RIA_SIDE_SIZE / 2: # if not filled but present ria must move top and gradient greater than atom count (side)
                to_choose.append((abs(atom.gradient.y), top_coord))
        if len(to_choose) == 0: # no swap required
            return None
        if len(to_choose) == 1: # either the top or bottow ria is invalid or only top/bottom requires to be swaped
            return to_choose[0]
        if to_choose[0][0] > to_choose[1][0]:
            return to_choose[0]
        else:
            return to_choose[1]

    def thermoswap_on_atom(self, atom: Atom):
        #print("thermo swap on atom")
        #print(atom.pos.x, atom.pos.y, atom.gradient.x, atom.gradient.y)
        hbswap = self.get_horizontal_best_swap(atom)
        vbswap = self.get_vertical_best_swap(atom)
        #print(hbswap, vbswap)
        if hbswap == None and vbswap == None: # no swap required
            return
        if hbswap == None and vbswap != None:
            self.switch_atoms(atom.pos.x, atom.pos.y, vbswap[1].x, vbswap[1].y)
            return
        if hbswap != None and vbswap == None:
            self.switch_atoms(atom.pos.x, atom.pos.y, hbswap[1].x, hbswap[1].y)
            return
        #print(hbswap[1].x, hbswap[1].y, vbswap[1].x, vbswap[1].y)
        if hbswap != None and vbswap != None and hbswap[0] > vbswap[0]:
            self.switch_atoms(atom.pos.x, atom.pos.y, hbswap[1].x, hbswap[1].y)
        elif hbswap != None and vbswap != None and hbswap[0] <= vbswap[0]:
            self.switch_atoms(atom.pos.x, atom.pos.y, vbswap[1].x, vbswap[1].y)
        else:
            raise Exception("unexpected bug")

    def thermoswap(self):
        for atom in self.atoms:
            self.thermoswap_on_atom(atom)



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
        target_a = Vec2(res_a[1].pos.x + res_a[1].gradient.x, res_a[1].pos.y + res_a[1].gradient.y)
        if res_b[0] == True:
            target_b = Vec2(res_b[1].pos.x + res_b[1].gradient.x, res_b[1].pos.y + res_b[1].gradient.y)
            res_a[1].pos.x, res_b[1].pos.x = res_b[1].pos.x, res_a[1].pos.x
            res_a[1].pos.y, res_b[1].pos.y = res_b[1].pos.y, res_a[1].pos.y
            res_b[1].gradient.set_xy(target_b.x - res_b[1].pos.x, target_b.y - res_b[1].pos.y)
        else:
            res_a[1].pos.x = x_b
            res_a[1].pos.y = y_b
        res_a[1].gradient.set_xy(target_a.x - res_a[1].pos.x, target_a.y - res_a[1].pos.y)
        self.rias[x_a][y_a], self.rias[x_b][y_b] = res_b[1], res_a[1]

    def apply_slight_shuffle_on_grid(self):
        ria_in_use: Ria = random.choice(self.ria_collected)
        new_pos = rand2Dexcept(ria_in_use.pos, GRID_SLIGHT_SHUFFLE_DISTANCE, self.len() - 1, 0, self.lenlen() - 1, 0)
        self.switch_rias(ria_in_use.pos.x, ria_in_use.pos.y, new_pos.x, new_pos.y)

    def apply_slight_shuffle_on_random_ria(self):
        ria_in_use: Ria = random.choice(self.ria_collected)
        ria_in_use.apply_slight_shuffle()

    def update_atom_gradients(self, previous, next):
        for ria in self.ria_collected:
            for atom in ria.atoms:
                atom.gradient.set_xy(0.0, 0.0)
                div = 0.0
                if previous != None:
                    div += 1.0
                    atom_previous_tuple = previous.find_atom(atom.id)
                    atom.gradient.x += (atom_previous_tuple[0].pos.x * RIA_SIDE_SIZE + atom_previous_tuple[1].pos.x) - (ria.pos.x * RIA_SIDE_SIZE + atom.pos.x)
                    atom.gradient.y += (atom_previous_tuple[0].pos.y * RIA_SIDE_SIZE + atom_previous_tuple[1].pos.y) - (ria.pos.y * RIA_SIDE_SIZE + atom.pos.y)
                if next != None:
                    div += 1.0
                    atom_next_tuple = next.find_atom(atom.id)
                    atom.gradient.x += (atom_next_tuple[0].pos.x * RIA_SIDE_SIZE + atom_next_tuple[1].pos.x) - (ria.pos.x * RIA_SIDE_SIZE + atom.pos.x)
                    atom.gradient.y += (atom_next_tuple[0].pos.y * RIA_SIDE_SIZE + atom_next_tuple[1].pos.y) - (ria.pos.y * RIA_SIDE_SIZE + atom.pos.y)
                atom.gradient.x /= div
                atom.gradient.y /= div

    def update_ria_gradients(self):
        for ria in self.ria_collected:
            ria.gradient.set_xy(0.0, 0.0)
            for atom in ria.atoms:
                ria.gradient.x += float(atom.gradient.x)
                ria.gradient.y += float(atom.gradient.y)
            ria.gradient.x /= float(len(ria.atoms))
            ria.gradient.y /= float(len(ria.atoms))

    def get_ria_from_pos(self, x: int, y: int): # "pos" in ria coords! (not absolute pos)
        if x < 0 or y < 0 or x >= self.len() or y >= self.lenlen(): # out of boundaries
            return (False, None)
        return (True, self.rias[x][y])

    def get_horizontal_best_swap(self, ria: Ria):
        left_coord = Vec2(ria.pos.x - 1, ria.pos.y)
        right_coord = Vec2(ria.pos.x + 1, ria.pos.y)
        left_tuple = self.get_ria_from_pos(left_coord.x, left_coord.y)
        right_tuple = self.get_ria_from_pos(right_coord.x, right_coord.y)
        to_choose: List[tuple[float, Vec2]] = []
        if left_tuple[0] == True: # left ria present
            if left_tuple[1] != None and left_tuple[1].gradient.x > ria.gradient.x: # if left ria filled and not ordered
                to_choose.append((abs(left_tuple[1].gradient.x - ria.gradient.x), left_coord))
            elif ria.gradient.x < 0 and abs(ria.gradient.x) > RIA_SIDE_SIZE / 2: # if not filled but present ria must move left
                to_choose.append((abs(ria.gradient.x), left_coord))
        if right_tuple[0] == True: # right ria present
            if right_tuple[1] != None and right_tuple[1].gradient.x < ria.gradient.x: # if right ria filled and not ordered
                to_choose.append((abs(right_tuple[1].gradient.x - ria.gradient.x), right_coord))
            elif ria.gradient.x > 0 and abs(ria.gradient.x) > RIA_SIDE_SIZE / 2: # if not filled but present ria must move right
                to_choose.append((abs(ria.gradient.x), right_coord))
        if len(to_choose) == 0: # no swap required
            return None
        if len(to_choose) == 1: # either the top or bottow ria is invalid or only top/bottom requires to be swaped
            return to_choose[0]
        if to_choose[0][0] > to_choose[1][0]:
            return to_choose[0]
        else:
            return to_choose[1]

    def get_vertical_best_swap(self, ria: Ria):
        bottom_coord = Vec2(ria.pos.x, ria.pos.y - 1)
        top_coord = Vec2(ria.pos.x, ria.pos.y + 1)
        bottom_tuple = self.get_ria_from_pos(bottom_coord.x, bottom_coord.y)
        top_tuple = self.get_ria_from_pos(top_coord.x, top_coord.y)
        to_choose: List[tuple[float, Vec2]] = []
        if bottom_tuple[0] == True: # bottom ria present
            if bottom_tuple[1] != None and bottom_tuple[1].gradient.y > ria.gradient.y: # if bottom ria filled and not ordered and diff grad greater than
                to_choose.append((abs(bottom_tuple[1].gradient.y - ria.gradient.y), bottom_coord))
            elif ria.gradient.y < 0 and abs(ria.gradient.y) > RIA_SIDE_SIZE / 2: # if not filled but the present ria must move bottom and gradient greater than atom count (side)
                to_choose.append((abs(ria.gradient.y), bottom_coord))
        if top_tuple[0] == True: # top ria present
            if top_tuple[1] != None and top_tuple[1].gradient.y < ria.gradient.y: # if top ria filled and not ordered
                to_choose.append((abs(top_tuple[1].gradient.y - ria.gradient.y), top_coord))
            elif ria.gradient.y > 0 and abs(ria.gradient.y) > RIA_SIDE_SIZE / 2: # if not filled but present ria must move top and gradient greater than atom count (side)
                to_choose.append((abs(ria.gradient.y), top_coord))
        if len(to_choose) == 0: # no swap required
            return None
        if len(to_choose) == 1: # either the top or bottow ria is invalid or only top/bottom requires to be swaped
            return to_choose[0]
        if to_choose[0][0] > to_choose[1][0]:
            return to_choose[0]
        else:
            return to_choose[1]

    def thermoswap_on_ria(self, ria: Ria):
        #print("thermo swap on ria")
        #print(ria.pos.x, ria.pos.y, ria.gradient.x, ria.gradient.y)
        hbswap = self.get_horizontal_best_swap(ria)
        vbswap = self.get_vertical_best_swap(ria)
        #print(hbswap, vbswap)
        if hbswap == None and vbswap == None: # no swap required
            return
        if hbswap == None and vbswap != None:
            self.switch_rias(ria.pos.x, ria.pos.y, vbswap[1].x, vbswap[1].y)
            return
        if hbswap != None and vbswap == None:
            self.switch_rias(ria.pos.x, ria.pos.y, hbswap[1].x, hbswap[1].y)
            return
        #print(hbswap[1].x, hbswap[1].y, vbswap[1].x, vbswap[1].y)
        if hbswap != None and vbswap != None and hbswap[0] > vbswap[0]:
            self.switch_rias(ria.pos.x, ria.pos.y, hbswap[1].x, hbswap[1].y)
        elif hbswap != None and vbswap != None and hbswap[0] <= vbswap[0]:
            self.switch_rias(ria.pos.x, ria.pos.y, vbswap[1].x, vbswap[1].y)
        else:
            raise Exception("unexpected bug")

    def apply_thermoswap_from_other_grid(self, previous, next):
        self.update_atom_gradients(previous, next)
        self.update_ria_gradients()
        for ria in self.ria_collected:
            self.thermoswap_on_ria(ria)
        #self.update_atom_gradients(previous, next)
        #for ria in self.ria_collected:
        #    ria.thermoswap()

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
        for stage in stages:
            self.family.append(generate_grid(ria_count, stage))

    def copy(self) -> Self:
        cp = copy.deepcopy(self)
        return cp

    def estimate_total_transfer_cost(self) -> int:
        total_cost = 0
        for i in range(len(self.family) - 1):
            total_cost += estimate_transfer_cost(self.family[i], self.family[i + 1])
        return total_cost

    def apply_thermoswap(self):
        for i in range(len(self.family)):
            self.family[i].apply_thermoswap_from_other_grid(self.family[i - 1] if i > 0 else None, self.family[i + 1] if i < len(self.family) - 1 else None)

    def shuffle(self, ria_pass_count: int, atom_pass_count: int):
        for grid in self.family:
            for i in range(ria_pass_count):
                grid.apply_slight_shuffle_on_grid()
            for j in range(atom_pass_count):
                grid.apply_slight_shuffle_on_random_ria()

    def print(self):
        for i in range(len(self.family)):
            print("Stage ", i)
            self.family[i].print()



class Group:

    def __init__(self, id, total_atom_count, current_atom_count) -> None:
        self.id = id
        self.total_atom_count = total_atom_count
        self.current_atom_count = current_atom_count



def generate_random_stage_affiliations(atom_count, max_atom_per_group):
    assert atom_count > max_atom_per_group and max_atom_per_group > 0
    groups: List[Group] = []
    total_atom_count = 0
    id = 0;
    while total_atom_count < atom_count:
        current_atom_count_for_group = random.choice(range(1, max_atom_per_group + 1 if (atom_count - total_atom_count) >= max_atom_per_group else (atom_count - total_atom_count + 1)))
        groups.append(Group(id, current_atom_count_for_group, 0));
        total_atom_count += current_atom_count_for_group;
        id += 1;
    atom_affiliations: List[int] = [0] * atom_count
    group_count = len(groups);
    for i in range(atom_count):
        group_index = random.choice(range(0, len(groups)))#crand_generate(0, (utype)groups.size() - 1);
        atom_affiliations[i] = groups[group_index].id;
        groups[group_index].current_atom_count += 1;
        if groups[group_index].current_atom_count >= groups[group_index].total_atom_count:
            del groups[group_index]
    return (group_count, atom_affiliations);



def convert_atom_affiliations_to_stage(group_count, atom_affiliations):
    stage: List[List[int]] = []
    for i in range(group_count):
        stage.append([])
    for i in range(len(atom_affiliations)):
        stage[atom_affiliations[i]].append(i)
    return stage



### TEST



RIA_COUNT = 50
#stage_0: List[List[int]] = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15], [16, 17], [18, 19], [20, 21], [22, 23], [24], [25], [26, 27, 28], [29, 30, 31], [32, 33, 34], [35]]
#stage_1: List[List[int]] = [[0, 1, 4, 5], [2, 3, 6, 7], [8, 9, 12, 13], [10, 11, 14, 15], [16, 20, 24], [17, 18, 19], [21, 22], [23], [25], [26, 29, 32, 35], [27, 30, 33], [28, 31, 34]]
#stage_2: List[List[int]] = [[0, 25, 8, 9], [1, 4, 6, 13], [2, 3, 10, 11], [5, 7, 12, 14], [15, 16], [20, 22], [17, 23], [24], [18], [19], [21], [27, 32, 33, 35], [26, 30, 34], [28, 29, 31]]

#RIA_COUNT = 4
#stage_0: List[List[int]] = [[0, 1, 2, 3], [4, 5, 6, 7]]
#stage_1: List[List[int]] = [[0, 4, 2, 6], [1, 5, 3, 7]]

stage_affiliation_0 = generate_random_stage_affiliations(1000, 4)
stage_affiliation_1 = generate_random_stage_affiliations(1000, 4)
stage_affiliation_2 = generate_random_stage_affiliations(1000, 4)
stage_0 = convert_atom_affiliations_to_stage(stage_affiliation_0[0], stage_affiliation_0[1])
stage_1 = convert_atom_affiliations_to_stage(stage_affiliation_1[0], stage_affiliation_1[1])
stage_2 = convert_atom_affiliations_to_stage(stage_affiliation_2[0], stage_affiliation_2[1])
print(stage_0)
print(stage_1)
print(stage_2)

family = Family(RIA_COUNT, [stage_0, stage_1, stage_2])
family.shuffle(20, 20)
family.print()
for i in range(200):
    print("Thermo swap: ", i, ", cost: ", family.estimate_total_transfer_cost())
    family.apply_thermoswap()
family.print()
