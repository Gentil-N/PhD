from enum import Enum
from typing import List
from typing import Self
from typing import Dict
import copy
import numpy as np
import random



def find_index_list(int_target: int, int_list_list: List[List[int]]) -> int:
    for i in range(len(int_list_list)):
        if int_target in int_list_list[i]:
            return i
    return -1



class Qbit:

    def __init__(self, id: int  = -1) -> None:
        self.id: int = id

    def set_id(self, id: int) -> None:
        if self.id != -1:
            raise Exception("qbit's id is already set")
        self.id = id



class Grid:

    def __init__(self, side_size: int) -> None:
        self.side_size = side_size
        dummy_qbit = Qbit()
        self.grid_array = [[dummy_qbit for i in range(self.side_size)] for j in range(self.side_size)]

    def find_qbit(self, id: int) -> tuple[int, int]:
        for i in range(self.side_size):
            for j in range(self.side_size):
                if self.grid_array[i][j] != None and self.grid_array[i][j].id == id:
                    return (i, j)
        return (-1, -1)

    def place_qbit(self, vpos: int, hpos: int, qbit: Qbit) -> None:
        if hpos < 0 or hpos >= self.side_size or vpos < 0 or vpos >= self.side_size:
            raise Exception("invalid location (" + str(hpos) + "," + str(vpos) + ") for grid(" + str (self.side_size) + ")")
        if self.grid_array[vpos][hpos] != None and self.grid_array[vpos][hpos].id != -1:
            raise Exception("trying to overwrite the existing qbit " + str(self.grid_array[vpos][hpos].id) + " with " + str(qbit.id))
        self.grid_array[vpos][hpos] = qbit

    def print_to_console(self) -> None:
        str_result = "Grid(" + str(self.side_size) + ")\n"
        for i in range(self.side_size):
            for j in range(self.side_size):
                qbit = self.grid_array[i][j]
                if qbit == None:
                    str_result += ".\t"
                elif qbit.id == -1:
                    str_result += ".\t"
                else:
                    str_result += str(qbit.id) + "\t"
            str_result += "\n"
        print(str_result)



class MorphQbit:

    def __init__(self, qbit = Qbit()) -> None:
        self.qbit = qbit
        self.gradient: tuple[int, int] = (0, 0)
    def set_qubit(self, qbit: Qbit) -> None:
        if self.qbit.id != -1:
            raise Exception("qbit inside morph qbit is already set")
        self.qbit = qbit



class MorphGrid:

    def __init__(self, grid_side_size: int) -> None:
        self.grid_side_size = grid_side_size
        self.morph_grid_side_size = self.grid_side_size * 2 + 1
        dummy_morph_qbit = MorphQbit()
        self.grid_array = [[dummy_morph_qbit for i in range(self.morph_grid_side_size)] for j in range(self.morph_grid_side_size)]

    def find_morph_qbit(self, id: int) -> tuple[MorphQbit, tuple[int, int]]:
        for i in range(self.morph_grid_side_size):
            for j in range(self.morph_grid_side_size):
                if self.grid_array[i][j] != None and self.grid_array[i][j].qbit.id == id:
                    return (self.grid_array[i][j], (i, j))
        return (MorphQbit(), (-1, -1))

    def _check_pos(self, pos: tuple[int, int]) -> bool:
        return not (pos[0] < 0 or pos[0] > len(self.grid_array) or pos[1] < 0 or pos[1] > len(self.grid_array[0]))

    def is_qbit_present(self, pos: tuple[int, int]) -> bool:
        if not self._check_pos(pos):
            raise Exception("invalid position (", pos[0], ",", pos[1], ")")
        if self.grid_array[pos[0]][pos[1]].qbit.id != -1:
            return True
        return False

    def _set_qbit_if_invalid(self, morph_qbit: MorphQbit, pos: tuple[int, int], shift: tuple[int, int]) -> bool:
        new_pos = (shift[0] + pos[0], shift[1] + pos[1])
        if self.is_qbit_present(new_pos):
            return False
        morph_qbit.gradient = (morph_qbit.gradient[0] - shift[0], morph_qbit.gradient[1] - shift[1])
        self.grid_array[new_pos[0]][new_pos[1]] = morph_qbit
        self.grid_array[pos[0]][pos[1]] = MorphQbit()
        return True

    def shift_qbit(self, pos: tuple[int, int] = (0, 0)) -> tuple[bool, tuple[int, int]]:
        if self.grid_array[pos[0]][pos[1]].qbit.id == -1:
            return (False, (0, 0))
        current_morph_qbit = self.grid_array[pos[0]][pos[1]]
        #shift = (-1, 0)
        #if self._set_qbit_if_invalid(current_morph_qbit, pos, shift):
        #    return (True, (pos[0] + shift[0], pos[1] + shift[1]))
        shift = (1, 0)
        if self._set_qbit_if_invalid(current_morph_qbit, pos, shift):
            return (True, (pos[0] + shift[0], pos[1] + shift[1]))
        #shift = (0, -1)
        #if self._set_qbit_if_invalid(current_morph_qbit, pos, shift):
        #    return (True, (pos[0] + shift[0], pos[1] + shift[1]))
        #shift = (0, 1)
        #if self._set_qbit_if_invalid(current_morph_qbit, pos, shift):
        #    return (True, (pos[0] + shift[0], pos[1] + shift[1]))
        raise Exception("cannot shift qbit ", current_morph_qbit.qbit.id, " at ", pos[0], " ", pos[1])

    def shift_all_targeted_qbit(self):
        for i in range(len(self.grid_array)):
            for j in range(len(self.grid_array[i])):
                current_morph_qbit = self.grid_array[i][j]
                targeted_qbit = (i + current_morph_qbit.gradient[0], j + current_morph_qbit.gradient[1])
                if targeted_qbit[0] == i and targeted_qbit[1] == j:
                    continue
                if self.grid_array[targeted_qbit[0]][targeted_qbit[1]].qbit.id == -1:
                    continue
                print("From ", current_morph_qbit.qbit.id, "(", i, ",", j, ") shifting ", self.grid_array[targeted_qbit[0]][targeted_qbit[1]].qbit.id, "(", targeted_qbit[0], ",", targeted_qbit[1], ") ")
                self.shift_qbit(targeted_qbit)
                #self.print_to_console()

    def _are_morph_qbit_compatible(self, qbit_a: MorphQbit, pos_a: tuple[int, int], qbit_b: MorphQbit, pos_b: tuple[int, int]) -> bool:
        new_pos_a: tuple[int, int] = (pos_a[0] + qbit_a.gradient[0], pos_a[1] + qbit_a.gradient[1])
        new_pos_b: tuple[int, int] = (pos_b[0] + qbit_b.gradient[0], pos_b[1] + qbit_b.gradient[1])
        if pos_a[0] < pos_b[0] and new_pos_a[0] < new_pos_b[0]:
            if pos_a[1] < pos_b[1] and new_pos_a[1] < new_pos_b[1]:
                return True
            if pos_a[1] > pos_b[1] and new_pos_a[1] > new_pos_b[1]:
                return True
            if pos_a[1] == pos_b[1] and new_pos_a[1] == new_pos_b[1]:
                return True
        if pos_a[0] > pos_b[0] and new_pos_a[0] > new_pos_b[0]:
            if pos_a[1] < pos_b[1] and new_pos_a[1] < new_pos_b[1]:
                return True
            if pos_a[1] > pos_b[1] and new_pos_a[1] > new_pos_b[1]:
                return True
            if pos_a[1] == pos_b[1] and new_pos_a[1] == new_pos_b[1]:
                return True
        if pos_a[0] == pos_b[0] and new_pos_a[0] == new_pos_b[0]:
            if pos_a[1] < pos_b[1] and new_pos_a[1] < new_pos_b[1]:
                return True
            if pos_a[1] > pos_b[1] and new_pos_a[1] > new_pos_b[1]:
                return True
            if pos_a[1] == pos_b[1] and new_pos_a[1] == new_pos_b[1]:
                raise Exception("unexpected bug: both qbits at same place? ", qbit_a.qbit.id, ", ", qbit_b.qbit.id)
        return False

    def _process_morph_qbit_into_groups(self, morph_qbit: MorphQbit, pos: tuple[int, int], compatibility_groups: List[List[tuple[tuple[int, int], MorphQbit]]], groups: List[List[int]]):
        for i in range(len(compatibility_groups)):
            compatibility_flag = True
            for elem in compatibility_groups[i]:
                if not self._are_morph_qbit_compatible(morph_qbit, pos, elem[1], elem[0]):
                    compatibility_flag = False
                    break
            if compatibility_flag:
                compatibility_groups[i].append((pos, morph_qbit))
                groups[i].append(morph_qbit.qbit.id)
                return # the qbit is placed in a group: the job is finished
        compatibility_groups.append([(pos, morph_qbit)]) # if no suitable group found: append a new group with the specific morph qbit
        groups.append([morph_qbit.qbit.id])

    def get_compatible_morph_qbit_move_groups(self) -> List[List[int]]:
        groups: List[List[int]] = []
        compatibility_groups: List[List[tuple[tuple[int, int], MorphQbit]]] = []
        for i in range(len(self.grid_array)):
            for j in range(len(self.grid_array[i])):
                if self.grid_array[i][j].qbit.id == -1 or (self.grid_array[i][j].gradient[0] == 0 and self.grid_array[i][j].gradient[1] == 0):
                    continue
                self._process_morph_qbit_into_groups(self.grid_array[i][j], (i, j), compatibility_groups, groups)
        return groups

    def select_next_move_group(self, groups: List[List[int]]) -> int:
        if len(groups) == 0:
            raise Exception("cannot select a specific group from an empty list")
        weights: List[int] = [0] * len(groups)
        for i in range(len(self.grid_array)):
            for j in range(len(self.grid_array)):
                if i % 2 == 1 or j % 2 == 0: # the qbits awaiting for morphing are necessarily on odd columnns and even lines
                    continue
                if self.grid_array[i - 1][j].qbit.id == -1: # if the qbit place above is empty, do not count the current one as a candidate to move
                    continue
                qbit_id: int = self.grid_array[i][j].qbit.id
                if qbit_id == -1:
                    continue
                index: int = find_index_list(qbit_id, groups)
                if index == -1:
                    raise Exception("unexcpected bug: qbit n°", qbit_id, " is not present in any moving group")
                weights[index] += 1
        max_index: int = weights.index(max(weights))
        #print("Weight list: ", len(weights))
        print("Weight list: ", weights)
        if weights[max_index] == 0:
            return 0
        return max_index

    def activate_move(self, id_list: List[int]) -> List:
        start_infos: List[tuple[MorphQbit, tuple[int, int]]] = []
        for id in id_list:
            start_infos.append(self.find_morph_qbit(id))
        agent_list: List[Agent] = []
        for i in range(len(start_infos)):
            vtarget = start_infos[i][1][0] + start_infos[i][0].gradient[0]
            htarget = start_infos[i][1][1] + start_infos[i][0].gradient[1]
            idtarget = self.grid_array[vtarget][htarget].qbit.id
            if idtarget != -1:
                raise Exception("position (", vtarget, ",", htarget, ") is busy with ", idtarget)
            agent_list.append(Agent(id_list[i], (start_infos[i][1][0], start_infos[i][1][1]), (vtarget, htarget)))
            #print("Agent: ", agent_list[-1].target, " ", (agent_list[-1].open_list.get(next(iter(agent_list[-1].open_list))) or Node((-1, -1), (-1, -1))).pos)
        # save and remove the current agents from the grid
        saved_morph_qbits = []
        for start_info in start_infos:
            saved_morph_qbits.append(self.grid_array[start_info[1][0]][start_info[1][1]])
            self.grid_array[start_info[1][0]][start_info[1][1]] = MorphQbit()
        # run the path finding for each agent
        reserve_table: ReserveTable = ReserveTable()
        for agent in agent_list:
            #print("Path finding:", agent.qbit_id, agent.init_pos, agent.target)
            agent.run(self, reserve_table)
        # replace the agents
        for i in range(len(start_infos)):
            self.grid_array[start_infos[i][1][0]][start_infos[i][1][1]] = saved_morph_qbits[i]
        #print("finished")
        return agent_list

    def print_to_console(self, show_gradient: bool = False) -> None:
        str_result = "MorphGrid(" + str(self.morph_grid_side_size) + ")\n"
        for i in range(self.morph_grid_side_size):
            for j in range(self.morph_grid_side_size):
                morph_qbit = self.grid_array[i][j]
                if morph_qbit == None:
                    str_result += "."
                elif morph_qbit.qbit.id == -1:
                    str_result += "."
                else:
                    str_result += str(morph_qbit.qbit.id)
                    if show_gradient:
                        str_result += "(" + str(morph_qbit.gradient[0]) + "," + str(morph_qbit.gradient[1]) + ")"
                str_result += "\t"
            str_result += "\n"
        print(str_result)



class Node:

    def __init__(self, pos: tuple[int, int], target: tuple[int, int], parent = None) -> None:
        self.parent = parent
        self.time_step = 0
        if self.parent != None:
            self.time_step = (parent or Node((-1, -1), (-1, -1))).time_step + 1
        self.pos: tuple[int, int] = pos
        self.g_score: int = 0
        self.h_score: int = 0
        self.f_score: int = 0
        self._compute_g_score()
        self._compute_h_score(target)
        self._compute_f_score()

    def _compute_g_score(self):
        if self.parent == None:
            self.g_score = 0
            return
        self.g_score = self.parent.g_score + abs(self.parent.pos[0] - self.pos[0]) + abs(self.parent.pos[1] - self.pos[1])

    def _compute_h_score(self, target: tuple[int, int]):
        self.h_score = abs(target[0] - self.pos[0]) + abs(target[1] - self.pos[1])

    def _compute_f_score(self):
        self.f_score = self.g_score + self.h_score



class ReserveTable:

    def __init__(self) -> None:
        self.time_steps: List[Dict[tuple[int, int], int]] = [] # dict(pos, qbit_id)
        self.last_targets: List[tuple[tuple[int, int], int]] = [] # tuple(pos, time_step)

    def add_path(self, path: List[Node], qbit_id: int):
        for ri in reversed(range(len(path))):
            time_step = len(path) - ri - 1
            if ri == 0: # then, this is the target reached by the qbit
                self.last_targets.append((path[ri].pos, time_step))
            if time_step == len(self.time_steps):
                self.time_steps.append({})
            if path[ri].pos in self.time_steps[time_step].keys():
                print("Before bug: ", ri, path[ri].pos, time_step, self.time_steps[time_step][path[ri].pos])
                raise Exception("unexpected bug: a position for a new path is already present in the reservation table")
            if time_step != path[ri].time_step:
                raise Exception("unexepcted bug: different time steps")
            self.time_steps[time_step][path[ri].pos] = qbit_id

    def is_position_colliding(self, current_pos: tuple[int, int], next_pos: tuple[int, int], next_time_step) -> bool:
        # Firstly, test if this is a target position already reached
        for pos_time_step_tuple in self.last_targets:
            if next_time_step > pos_time_step_tuple[1] and next_pos[0] == pos_time_step_tuple[0][0] and next_pos[1] == pos_time_step_tuple[0][1]:
                return True
        # Secondly, test if the position is reserved by another travelling qbit
        if next_time_step < len(self.time_steps):
            # Vertex collision test
            next_dict = self.time_steps[next_time_step]
            if next_pos in next_dict.keys():
                return True
            # Edge collision test
            if next_time_step > 0:
                current_dict: Dict[tuple[int, int], int] = self.time_steps[next_time_step - 1]
                if next_pos in current_dict.keys(): # then, next pos already occupied at previous time step
                    qbit_id: int = current_dict[next_pos]
                    if current_pos in next_dict.keys() and next_dict[current_pos] == qbit_id: # then, qbits are swaped
                        return True
        return False

class Agent:

    def __init__(self, qbit_id: int, init_pos: tuple[int, int], target: tuple[int, int]) -> None:
        self.qbit_id = qbit_id
        self.init_pos = init_pos
        self.closed_list: Dict[tuple[int, int], Node] = {}
        self.open_list: Dict[tuple[int, int], Node] = {init_pos: Node(init_pos, target)}
        self.target: tuple[int, int] = target
        self.target_reached: bool = False
        self.path: List[Node] = []

    def _build_path(self):
        if not self.target_reached:
            raise Exception("tempting to build path whereas the target is not reached")
        if not(self.target in self.closed_list.keys()):
            raise Exception("target is reached but not evaluated in the closed list")
        current_node = self.closed_list.get(self.target) or Node((-1, -1), (-1, -1))
        while current_node != None:
            self.path.append(current_node)
            current_node = current_node.parent

    def _check_push_node(self, next_pos: tuple[int, int], parent_node: Node, morph_grid: MorphGrid, reserve_table: ReserveTable):
        #print("checking ", new_pos)
        if next_pos[0] < 0 or next_pos[0] >= morph_grid.morph_grid_side_size or next_pos[1] < 0 or next_pos[1] >= morph_grid.morph_grid_side_size: # outside of the grid
            return
        if morph_grid.is_qbit_present(next_pos): # this is an obstacle: cannot add it to the open list
            return
        if next_pos in self.closed_list.keys() or next_pos in self.open_list.keys(): # already in open or closed list
            return
        if reserve_table.is_position_colliding(parent_node.pos, next_pos, parent_node.time_step + 1): # colliding with other travelling qbit
            return
        #print(next_pos, " added!")
        new_node = Node(pos=next_pos, target=self.target, parent=parent_node)
        self.open_list[new_node.pos] = new_node

    def run(self, morph_grid: MorphGrid, reserve_table: ReserveTable):
        while not self.target_reached:
            if len(self.open_list) == 0:
                raise Exception("failed to reach the target")
            if self.target in self.open_list.keys():
                #print("found!")
                self.closed_list[self.target] = self.open_list[self.target]
                self.open_list.pop(self.target)
                self.target_reached = True
                self._build_path()
                reserve_table.add_path(self.path, self.qbit_id)
                return
            min_f_score_node: Node = self.open_list.get(next(iter(self.open_list))) or Node((-1, -1), (-1, -1))
            #print("Agent: ", self.target, " min pos ", min_f_score_node.pos, " f_score ", min_f_score_node.f_score)
            if min_f_score_node.pos[0] == -1 or min_f_score_node.pos[1] == -1:
                raise Exception("unable to get the first node of the open list")
            for item in self.open_list.items():
                #print("item ", item[1].pos, " f_score ", item[1].f_score)
                if min_f_score_node.f_score > item[1].f_score:
                    min_f_score_node = item[1]
            #print("min chosen ", min_f_score_node.pos, " f_score ", min_f_score_node.f_score)
            self.closed_list[min_f_score_node.pos] = min_f_score_node
            self.open_list.pop(min_f_score_node.pos)
            self._check_push_node((min_f_score_node.pos[0] - 1, min_f_score_node.pos[1]), min_f_score_node, morph_grid, reserve_table)
            self._check_push_node((min_f_score_node.pos[0] + 1, min_f_score_node.pos[1]), min_f_score_node, morph_grid, reserve_table)
            self._check_push_node((min_f_score_node.pos[0], min_f_score_node.pos[1] + 1), min_f_score_node, morph_grid, reserve_table)
            self._check_push_node((min_f_score_node.pos[0], min_f_score_node.pos[1] - 1), min_f_score_node, morph_grid, reserve_table)
            #for i in self.open_list.items():
            #    print(i[1].pos, i[1].f_score)



def create_morph_grid(init_grid: Grid, target_grid: Grid) -> MorphGrid:
    if init_grid.side_size != target_grid.side_size:
        raise Exception("both grids must be identicals")
    grid_side_size = init_grid.side_size
    morph_grid = MorphGrid(grid_side_size)
    for i in range(init_grid.side_size):
        for j in range(init_grid.side_size):
            if init_grid.grid_array[i][j] != None and init_grid.grid_array[i][j].id != -1:
                qbit = init_grid.grid_array[i][j]
                target_pos = target_grid.find_qbit(init_grid.grid_array[i][j].id)
                if target_pos[0] == -1 or target_pos[1] == -1:
                    raise Exception("qbit " + str(qbit.id) + " is missing in the target grid")
                morph_qbit = MorphQbit(qbit)
                morph_qbit.gradient = ((target_pos[0] - i) * 2, (target_pos[1] - j) * 2)
                morph_grid.grid_array[i * 2 + 1][j * 2 + 1] = morph_qbit
    return morph_grid



def animate_move(previous_morph_grid: MorphGrid, agent_list: List[Agent]):
    morph_grid = copy.deepcopy(previous_morph_grid)
    step = 0
    change = True
    previous_pos: List[tuple[int, int]] = []
    for agent in agent_list:
        if not agent.target_reached:
            raise Exception("target not reached: cannot process animation")
        if len(agent.path) == 0:
            continue
        previous_pos.append(agent.init_pos)
    while change:
        change = False
        for i in range(len(agent_list)):
            #print("agent ", i)
            if not agent_list[i].target_reached:
                raise Exception("target not reached: cannot process animation")
            if len(agent_list[i].path) == 0:
                continue
            qbit: MorphQbit = morph_grid.grid_array[previous_pos[i][0]][previous_pos[i][1]]
            morph_grid.grid_array[previous_pos[i][0]][previous_pos[i][1]] = MorphQbit()
            #print(step, " ", len(agent_list[i].path))
            if step < len(agent_list[i].path):
                change = True
                node: Node = agent_list[i].path[len(agent_list[i].path) - 1 - step]
            else:
                node: Node = agent_list[i].path[0]
            qbit.gradient = (previous_pos[i][0] + qbit.gradient[0] - node.pos[0], previous_pos[i][1] + qbit.gradient[1] - node.pos[1])
            previous_pos[i] = node.pos
            morph_grid.grid_array[node.pos[0]][node.pos[1]] = qbit
        print("Step n°", step)
        morph_grid.print_to_console()
        step += 1
    return morph_grid



def build_init_target_grid(num_side_qbits: int) -> tuple[Grid, Grid]:
    num_qbits: int = num_side_qbits**2
    qbit_list: List[Qbit] = []
    for i in range(num_qbits):
        qbit_list.append(Qbit(i))
    init_grid: Grid = Grid(num_side_qbits)
    for i in range(num_side_qbits):
        for j in range(num_side_qbits):
            init_grid.place_qbit(i, j, qbit_list[i * num_side_qbits + j])
    target_grid: Grid = Grid(num_side_qbits)
    for i in range(num_side_qbits):
        for j in range(num_side_qbits):
            index = qbit_list.index(random.choice(qbit_list))
            target_grid.place_qbit(i, j, qbit_list[index])
            del qbit_list[index]
    return (init_grid, target_grid)



def check_collision(agent_list: List[Agent]) -> int:
    change = True
    index = 0
    num_collision = 0
    while change:
        for i in range(len(agent_list)):
            for j in range(i, len(agent_list)):
                if i == j:
                    continue
                index_i = len(agent_list[i].path) - 1 - min(index, len(agent_list[i].path) - 1)
                index_j = len(agent_list[j].path) - 1 - min(index, len(agent_list[j].path) - 1)
                #print(index_i, index_j, index, len(agent_list[i].path), len(agent_list[j].path))
                if agent_list[i].path[index_i].pos[0] ==  agent_list[j].path[index_j].pos[0] and agent_list[i].path[index_i].pos[1] ==  agent_list[j].path[index_j].pos[1]:
                    #print(i, j, index_i, index_j, agent_list[i].init_pos[0], agent_list[i].init_pos[1], agent_list[j].init_pos[0], agent_list[j].init_pos[1], agent_list[i].path[index_i].pos[0], agent_list[j].path[index_j].pos[0], agent_list[i].path[index_i].pos[1], agent_list[j].path[index_j].pos[1])
                    #print("agent ", i)
                    #for node in agent_list[i].path:
                    #    print(node.pos[0], node.pos[1])
                    #print("agent ", j)
                    #for node in agent_list[j].path:
                    #    print(node.pos[0], node.pos[1])
                    num_collision += 1
        index += 1
        change = False
        for agent in agent_list:
            if index < len(agent.path):
                change = True
                break
    return num_collision



def run_random_morphing(num_side_qbits: int) -> int:
    init_target_grids: tuple[Grid, Grid] = build_init_target_grid(num_side_qbits)
    init_target_grids[0].print_to_console()
    init_target_grids[1].print_to_console()
    morph_grid = create_morph_grid(init_target_grids[0], init_target_grids[1])
    print("Original morph grid")
    morph_grid.print_to_console()
    morph_grid.shift_all_targeted_qbit()
    print("Shifted morph grid")
    morph_grid.print_to_console()
    groups = morph_grid.get_compatible_morph_qbit_move_groups()
    print("Morphing starting")
    step_count = 0
    total_collision_count = 0
    while len(groups):
        i: int = morph_grid.select_next_move_group(groups)
        print("Group n°", i, " in process: ", groups[i])
        agent_list = morph_grid.activate_move(groups[i])
        group_collision_count = check_collision(agent_list)
        total_collision_count += group_collision_count
        print("Collision count: ", group_collision_count)
        morph_grid = animate_move(morph_grid, agent_list)
        del groups[i]
        step_count += 1
    print("Morphing terminated in ", step_count, " morph steps including ", total_collision_count, " collisions and averaging ", num_side_qbits**2 / step_count, "qbits per group")
    return total_collision_count


### TEST



total_collision_count = run_random_morphing(10)
