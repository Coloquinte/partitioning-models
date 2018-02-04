#!/usr/bin/python3

import localsolver
import sys

def read_graph(graph_file_name):
    edges = []
    edge_weights = []
    node_weights = []
    
    with open(graph_file_name) as f:
        n_edges, n_nodes, mode = [int (n) for n in f.readline().split()]
        # .hgr magic number
        if not mode in [0, 1, 10, 11]:
            raise Exception("Invalid hmetis mode")
        has_edge_weights = mode in [1, 11]
        has_node_weights = mode in [10, 11]
        # edges
        for i in range(n_edges):
            pins = [int(n) - 1 for n in f.readline().split()]
            weight = 1
            if has_edge_weights:
                weight = pins[0]
                pins = pins[1:]
            for p in pins:
                p = p-1
            edges.append(pins)
            edge_weights.append(weight)
        # node weights
        for i in range(n_nodes):
            if has_node_weights:
                node_weights.append(int(f.readline().split()[0]))
            else:
                node_weights.append(1)
    return (edges, edge_weights, node_weights)

def check_graph(edges, edge_weights, node_weights):
    for e in edges:
        for p in e:
            assert p >= 0 and p < len(node_weights)
    for w in edge_weights:
        assert w >= 0
    for w in node_weights:
        assert w >= 0

class ModelBuilder(object):
    def __init__(self, edges, edge_weights, node_weights):
        self.edges = edges
        self.edge_weights = edge_weights
        self.node_weights = node_weights
        self.ls = localsolver.LocalSolver()
        self.model = self.ls.get_model()
        self.n_parts = 2
        self.margin = 5.0
        self.node_placement = []
        self.edge_counters = []
        self.edge_degrees = []
        self.solution = None

    def build(self):
        self.build_node_variables()
        self.build_node_constraints()
        self.build_edge_counters()
        self.build_cut_cost()

    def solve(self, time_limit=None, iteration_limit=None):
        self.model.close()
        if iteration_limit == None:
            iteration_limit = 100*len(self.node_weights)
        self.ls.create_phase().iteration_limit = iteration_limit
        self.ls.solve()
        self.solution = []
        for i in range(len(self.node_placement)):
            pos = -1
            placement = self.node_placement[i]
            for j in range(self.n_parts):
                if placement[j].value:
                    pos = j
            self.solution.append(pos)

    def build_node_variables(self):
        self.node_placement = []
        for w in self.node_weights:
            places = [self.model.bool() for i in range(self.n_parts)]
            self.node_placement.append(places)
            self.model.add_constraint(self.model.sum(places) == 1)

    def build_node_constraints(self):
        for part in range(self.n_parts):
            self.build_node_constraint(part)

    def build_node_constraint(self, part):
        tot_weight = sum(self.node_weights)
        weight_per_part = int(tot_weight * (1.0 + self.margin / 100.0) / self.n_parts)
        weights_on_part = []
        for i in range(len(self.node_weights)):
            decision = self.node_weights[i]
            weight = self.node_placement[i][part]
            weights_on_part.append(decision * weight)
        self.model.add_constraint(self.model.sum(weights_on_part) < weight_per_part)

    def build_edge_counters(self):
        for i in range(len(self.edges)):
            pins = self.edges[i]
            counters = []
            for j in range(self.n_parts):
                counters.append(self.model.sum([self.node_placement[p][j] for p in pins]))
            self.edge_counters.append(counters)
            degree = self.model.sum(counters)
            self.edge_degrees.append(degree)

    def build_cut_cost(self):
        cut = self.model.sum([d >= 2 for d in self.edge_degrees])
        self.model.minimize(cut)


if len(sys.argv) < 2:
    print("Usage: ls_partitioning.py graph.hgr")
    sys.exit(1)

graph_file_name = sys.argv[1]

edges, edge_weights, node_weights = read_graph(graph_file_name)
check_graph(edges, edge_weights, node_weights)

builder = ModelBuilder(edges, edge_weights, node_weights)
builder.build()
builder.solve()

with localsolver.LocalSolver() as ls:
    print("Stuff");

