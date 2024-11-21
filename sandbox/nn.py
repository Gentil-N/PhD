import networkx as nx
import torch
from torch_geometric.nn import GCNConv
from torch_geometric.data import Data
from torch.utils.data import DataLoader, TensorDataset
#from torch_geometric.loader import DataLoader
from torch_geometric.datasets import TUDataset  # For graph data sets
from torch_geometric.utils import from_networkx
import torch.nn.functional as F
from torch.optim.lr_scheduler import StepLR, LambdaLR
from torch_geometric.nn import GCNConv
from typing import List
from typing import Self
from typing import Dict
import copy
import numpy as np
import random
import math



def read_log_file(file_name: str):
    file = open(file_name, "r")
    header_lines = []
    gene_lines = []
    current_lines = header_lines
    for line in file.readlines():
        current_lines.append(line)
        if "HEADER END" in line:
            current_lines = gene_lines
    file.close()

    header_dict = {}
    for line in header_lines:
        if "=" in line:
            splitted = line.split("=")
            val = splitted[1].strip()
            if "atomaffiliation" in splitted[0]:
                nums = val.split(",")
                val = []
                for num in nums:
                    val.append(int(num))
            header_dict[splitted[0].strip()] = val

    gene_dicts = []
    current_dict = None
    for line in gene_lines:
        if "GENE END" in line:
            gene_dicts.append(current_dict)
            current_dict = None
        if current_dict != None and "=" in line:
            splitted = line.split("=")
            val = splitted[1].strip()
            if "_" in splitted[0]:
                coord_cut = splitted[1].split(",")
                val = (int(coord_cut[0].strip()), int(coord_cut[1].strip()))
            current_dict[splitted[0].strip()] = val
        if "GENE START" in line:
            current_dict = {}

    return [header_dict, gene_dicts]



def create_connector(header_dict):
    atom_count = int(header_dict["atomcount"])
    stage_count = int(header_dict["stagecount"])
    connector = [] # one connector per stage
    for i in range(stage_count):
        atom_affiliations = header_dict["atomaffiliations" + str(i)]
        current_connector = []
        for j in range(len(atom_affiliations)):
            atom_group_id = atom_affiliations[j]
            for k in range(j + 1, len(atom_affiliations)):
                mask = 0.0
                if atom_group_id == atom_affiliations[k]:
                    mask = 1.0
                current_connector.append(mask)
        connector.extend(current_connector)
    return connector



def create_dissimilarities(header_dict, gene_dict):
    atom_count = int(header_dict["atomcount"])
    stage_count = int(header_dict["stagecount"])
    dissimilarities = [] # one dissimilarity (~matrix in list shape) per stage
    #for i in range(stage_count):
    #    #print("stage ", i)
    #    current_dissimilarity = []
    #    for j in range(atom_count):
    #        current_tuple_data = gene_dict[str(i) + "_" + str(j)]
    #        for k in range(j + 1, atom_count):
    #            #print(j, k)
    #            next_tuple_data = gene_dict[str(i) + "_" + str(k)]
    #            current_dissimilarity.append(
    #                math.sqrt(float(current_tuple_data[0] - next_tuple_data[0])**2 + float(current_tuple_data[1] - next_tuple_data[1])**2)
    #            )
    #    dissimilarities.append(current_dissimilarity);
    for i in range(stage_count):
        #print("stage ", i)
        current_dissimilarity = []
        for j in range(atom_count):
            row_data = []
            current_tuple_data = gene_dict[str(i) + "_" + str(j)]
            for k in range(atom_count):
                #print(j, k)
                next_tuple_data = gene_dict[str(i) + "_" + str(k)]
                row_data.append(
                    math.sqrt(float(current_tuple_data[0] - next_tuple_data[0])**2 + float(current_tuple_data[1] - next_tuple_data[1])**2)
                )
            current_dissimilarity.append(row_data)
        dissimilarities.append(current_dissimilarity)
    return dissimilarities



def load_nn_training_data(file):
    file_result = read_log_file(file)
    connector = create_connector(file_result[0])
    dissimilarities = create_dissimilarities(file_result[0], file_result[1][0])
    training_data_inputs = []
    training_data_outputs = []
    for i in range(len(dissimilarities)): # for each stage
        for j in range(len(dissimilarities[i])): # for each row inside the matrix
            for k in range(len(dissimilarities[i][j])): # for each value inside a given row
                dissimilarity_value = dissimilarities[i][j][k]
                in_data = [float(i), float(j), float(k)]
                in_data.extend(connector.copy())
                training_data_inputs.append(in_data)
                training_data_outputs.append([dissimilarity_value])
    return (training_data_inputs, training_data_outputs)



# Define the neural network
class SimpleNN(torch.nn.Module):
    def __init__(self, in_neuron_count, hidden_neuron_count_per_layer, out_neuron_count):
        super(SimpleNN, self).__init__()
        self.fc1 = torch.nn.Linear(in_neuron_count, hidden_neuron_count_per_layer)
        self.fc2 = torch.nn.Linear(hidden_neuron_count_per_layer, hidden_neuron_count_per_layer)
        self.fc3 = torch.nn.Linear(hidden_neuron_count_per_layer, hidden_neuron_count_per_layer)
        self.fc4 = torch.nn.Linear(hidden_neuron_count_per_layer, out_neuron_count)

    def forward(self, x):
        x = torch.sigmoid(self.fc1(x))
        x = torch.sigmoid(self.fc2(x))
        x = torch.sigmoid(self.fc3(x))
        x = self.fc4(x)
        return x



def create_graph(header_dict, gene_dict):
    graph_target = nx.Graph()
    atom_count = int(header_dict["atomcount"])
    stage_count = int(header_dict["stagecount"])
    for i in range(atom_count):
        one_hot_position_stage_default = []
        for j in range(stage_count):
            tuple_data = gene_dict[str(j) + "_" + str(i)]
            one_hot_position_stage_default.extend([float(tuple_data[0]), float(tuple_data[1])])
        graph_target.add_node(i)
        nx.set_node_attributes(graph_target, {i : {"one_hot_position_stage" : one_hot_position_stage_default}})
        #print(i, graph.nodes[i]["one_hot_position_stage"])

    stages = []
    for i in range(stage_count):
        current_stage = []
        atom_affiliations = header_dict["atomaffiliations" + str(i)]
        for _ in range(int(header_dict["groupcount" + str(i)])):
            current_stage.append([])
        for j in range(len(atom_affiliations)):
            current_stage[atom_affiliations[j]].append(j)
        stages.append(current_stage)
    for i in range(len(stages)):
        for group in stages[i]:
            for k in range(len(group) - 1):
                for j in range(k + 1, len(group)):
                    if not graph_target.has_edge(group[k], group[j]):
                        one_hot_active_stage_default = [0.0] * stage_count
                        graph_target.add_edge(group[k], group[j])
                        nx.set_edge_attributes(graph_target, {(group[k], group[j]) : {"one_hot_active_stage" : one_hot_active_stage_default}})
                    graph_target.edges[group[k], group[j]]["one_hot_active_stage"][i] = 1.0
                    #print(group[k], group[j], graph.edges[group[k], group[j]])
    return graph_target



def load_gnn_data_from_file(file):
    print("Loading ", file)
    result = read_log_file(file)
    graph = create_graph(result[0], result[1][0])
    data = from_networkx(graph, group_node_attrs=["one_hot_position_stage"], group_edge_attrs=["one_hot_active_stage"])
    #print(data.edge_index)
    #print(data.x)
    #print(data.edge_attr)
    return Data(x=torch.rand(data.x.size(), dtype=torch.float) * 0.1, edge_index=data.edge_index, edge_attr=data.edge_attr, y=data.x)
    #return Data(x=data.x, edge_index=data.edge_index, edge_attr=data.edge_attr, y=data.x)



# Example GCN model with two GCN layers
class GCN(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super(GCN, self).__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, hidden_channels)
        self.conv3 = GCNConv(hidden_channels, out_channels)
        #self.conv4 = GCNConv(hidden_channels, out_channels)
        #self.conv5 = GCNConv(hidden_channels, hidden_channels)
        #self.conv6 = GCNConv(hidden_channels, hidden_channels)
        #self.conv7 = GCNConv(hidden_channels, out_channels)

    def forward(self, x, edge_index, edge_attr):
        x = self.conv1(x, edge_index)
        x = F.sigmoid(x)
        x = self.conv2(x, edge_index)
        x = F.sigmoid(x)
        x = self.conv3(x, edge_index)
        #x = F.relu(x)
        #x = self.conv4(x, edge_index)
        #x = F.relu(x)
        #x = self.conv5(x, edge_index)
        #x = F.relu(x)
        #x = self.conv6(x, edge_index)
        #x = F.relu(x)
        #x = self.conv7(x, edge_index)
        return x



def create_fake_data():
    node_features = torch.tensor([
            [0.1, 0.2, 0.2, 0.1],  # Node 0 features
            [0.3, 0.4, 0.4, 0.5],  # Node 1 features
            [0.5, 0.6, 0.6, 0.5]   # Node 2 features
        ], dtype=torch.float)

    edge_index = torch.tensor([
        [0, 0, 1, 1],  # Source nodes
        [1, 1, 2, 2]   # Target nodes
    ], dtype=torch.long)

    edge_attr = torch.tensor([
        [1.0, 0.0],  # Edge feature
        [1.0, 1.0],
        [0.0, 1.0],
        [1.0, 1.0]
    ], dtype=torch.float)

    # Assuming output node positions for training
    target_node_positions = torch.tensor([
        [0.0, 0.1, 0.3, 0.5],
        [0.4, 0.5, 0.4, 0.6],
        [0.3, 0.4, 0.6, 0.7]
    ], dtype=torch.float)

    return Data(x=node_features, edge_index=edge_index, edge_attr=edge_attr, y=target_node_positions)



### TEST n°1



print("Init model")
training_data_inputs = []
training_data_outputs = []
for i in range(20):
    data_tuple = load_nn_training_data("./build/release/log" + str(i))
    training_data_inputs.extend(data_tuple[0])
    training_data_outputs.extend(data_tuple[1])

#print(len(training_data_inputs))
train_dataset = TensorDataset(torch.tensor(training_data_inputs, dtype=torch.float), torch.tensor(training_data_outputs, dtype=torch.float))
train_loader = DataLoader(train_dataset, batch_size=10, shuffle=True)
# Initialize the model, loss function, and optimizer
in_neuron_count = len(training_data_inputs[0])
model = SimpleNN(in_neuron_count, 2*in_neuron_count, 1)
criterion = torch.nn.MSELoss()  # Mean Squared Error Loss for regression
optimizer = torch.optim.SGD(model.parameters(), lr=0.01)

print("Start training")
# Training loop
num_epochs = 1000
total_loss = 10.0
epoch = 0
while total_loss > 0.07:
    total_loss = 0.0
    for batch_inputs, batch_targets in train_loader:
        optimizer.zero_grad()               # Reset gradients
        output = model(batch_inputs)                   # Forward pass
        loss = criterion(output, batch_targets)         # Compute loss
        loss.backward()                     # Backward pass
        optimizer.step()                    # Update weights
        total_loss += loss.item()
    total_loss = float(total_loss) / float(len(train_loader))
    print("epoch", epoch, "loss", total_loss)
    epoch += 1

with torch.no_grad():
    #for i in range(min(len(training_data_inputs), len(training_data_outputs))):
    #    test_output = model(torch.tensor(training_data_inputs[i]))
    #    print("exact value: ", training_data_outputs[i], "prediction: ", test_output.item())
    print("Genralization?")
    validation_data_inputs = []
    validation_data_output = []
    for i in range(3):
        data_tuple = load_nn_training_data("./build/release/log70" + str(i))
        validation_data_inputs.extend(data_tuple[0])
        validation_data_output.extend(data_tuple[1])
    for i in range(min(len(validation_data_inputs), len(validation_data_output))):
        test_output = model(torch.tensor(validation_data_inputs[i]))
        print("exact value: ", validation_data_output[i], "prediction: ", test_output.item())

exit()



### TEST n°2



#log_files = ["./log-test"]
#dataset = [load_gnn_data_from_file(log_files[0]) for _ in range(100)]
#print(dataset)
dataset = []
for i in range(100):
    dataset.append(load_gnn_data_from_file("./build/release/log" + str(i)))
#dataloader = DataLoader(dataset, batch_size=10, shuffle=True)
#dataloader = DataLoader(dataset, batch_size=1, shuffle=True) # conflict with other DataLoader
print("Data set fully loaded")

# Initialize the GCN model
model = GCN(in_channels=6, hidden_channels=80, out_channels=6)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
loss_fn = torch.nn.MSELoss()
#lambda_group2 = lambda epoch: 0.95 ** epoch
#scheduler = LambdaLR(optimizer, lr_lambda=lambda2)
#scheduler = StepLR(optimizer, step_size=100, gamma=0.00001)
print("GCN model initialized")

# Training loop
average_loss = 10
i = 0
while(average_loss > 1.0):  # Set the number of epochs
    total_loss = 0
    for batch in dataloader:
        if i > 20000:
            optimizer.param_groups[0]['lr'] = 0.005
        if i > 50000:
            optimizer.param_groups[0]['lr'] = 0.001
        if i > 100000:
            optimizer.param_groups[0]['lr'] = 0.0005
        optimizer.zero_grad()  # Clear gradients
        # Forward pass
        #print(batch.y)
        out = model(batch.x, batch.edge_index, batch.edge_attr)
        # Compute loss (MSE loss with respect to target positions)
        loss = loss_fn(out, batch.y)
        loss.backward()  # Backward pass
        optimizer.step()  # Update model parameters
        #scheduler.step(epoch=i)
        total_loss += loss.item() * batch.num_graphs  # Accumulate batch loss
        average_loss = total_loss / len(dataset)
    i += 1
    print(f"epoch {i} loss: {average_loss:.4f}")

# Test model on a new graph (or evaluate performance as needed)
test_data_set = [load_gnn_data_from_file("./build/release/log0"), load_gnn_data_from_file("./build/release/log10"), load_gnn_data_from_file("./build/release/log20")]
with torch.no_grad():
    for test_data in test_data_set:
        predicted_positions = model(test_data.x, test_data.edge_index, test_data.edge_attr)
        print("Predicted Node Positions:\n", predicted_positions)

torch.save(model.state_dict(), "./gnn-comp-model")
