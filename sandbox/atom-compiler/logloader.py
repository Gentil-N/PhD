import networkx as nx
from typing import List
from typing import Self
from typing import Dict
import copy
import numpy as np
import random
import math

def read_log_file(file_name: str):
    file = open(file_name, "r")
    header_dict = {}
    for line in file.readlines():
        if "=" in line:
            splited = line.split("=")
            header_dict[splited[0]] = splited[1]
    file.close()
