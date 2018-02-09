#!/usr/bin/python3

import localsolver
import sys
import random
import math

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
        self.node_placement = None
        self.edge_counters = None
        self.edge_degrees = None
        self.solution = None
        self.initial_solution = None

    def check(self):
        assert len(self.edges) == len(self.edge_weights)
        if self.node_placement != None:
            assert len(self.node_placement) == len(self.node_weights)
        if self.edge_counters != None:
            assert len(self.edge_counters) == len(self.edge_weights)
        if self.edge_degrees != None:
            assert len(self.edge_degrees) == len(self.edge_weights)
        if self.solution != None:
            assert len(self.solution) == len(self.node_weights)
        if self.initial_solution != None:
            assert len(self.initial_solution) == len(self.node_weights)

    def build(self):
        self.build_node_variables()
        self.build_node_constraints()
        self.build_edge_counters()
        self.build_cut_cost()

    def solve(self, time_limit=None, iteration_limit=None, seed=0):
        self.ls.get_param().set_verbosity(0)
        self.model.close()
        self.init_placement()
        if iteration_limit == None:
            #iteration_limit = 100*len(self.node_weights)
            iteration_limit = 1000000
        self.ls.create_phase().iteration_limit = iteration_limit
        self.ls.get_param().set_seed(seed)
        self.ls.solve()
        self.read_solution()

    def read_solution(self):
        self.solution = []
        for i in range(len(self.node_placement)):
            pos = -1
            placement = self.node_placement[i]
            for j in range(self.n_parts):
                if placement[j].value:
                    pos = j
            self.solution.append(pos)
        self.objective_value = self.model.get_objective(0).value

    def build_node_variables(self):
        self.node_placement = []
        for w in self.node_weights:
            places = [self.model.bool() for i in range(self.n_parts)]
            self.node_placement.append(places)
            self.model.add_constraint(self.model.sum(places) >= 1)

    def init_placement(self):
        if self.initial_solution is None:
            self.initial_solution = [random.randrange(self.n_parts) for i in self.node_placement]

        for i in range(len(self.node_placement)):
            places = self.node_placement[i]
            for p in places:
                p.set_value(False)
            init = self.initial_solution[i]
            places[init].set_value(True)

    def build_node_constraints(self):
        tot_weight = sum(self.node_weights)
        self.weight_per_part = int(tot_weight * (1.0 + self.margin / 100.0) / self.n_parts)
        for part in range(self.n_parts):
            self.build_node_constraint(part)

    def build_node_constraint(self, part):
        weights_on_part = []
        for i in range(len(self.node_weights)):
            decision = self.node_weights[i]
            weight = self.node_placement[i][part]
            weights_on_part.append(decision * weight)
        self.model.add_constraint(self.model.sum(weights_on_part) < self.weight_per_part)

    def build_edge_counters(self):
        self.edge_counters = []
        self.edge_degrees = []
        for i in range(len(self.edges)):
            pins = self.edges[i]
            counters = []
            for j in range(self.n_parts):
                counters.append(self.model.sum([self.node_placement[p][j] for p in pins]))
            self.edge_counters.append(counters)
            occupied = [c != 0 for c in counters]
            degree = self.model.sum(occupied)
            self.edge_degrees.append(degree)

    def build_cut_cost(self):
        edge_costs = [self.edge_weights[i] * (self.edge_degrees[i] >= 2) for i in range(len(self.edges))]
        cut = self.model.sum(edge_costs)
        self.model.minimize(cut)

class MultilevelSolver(object):
    def __init__(self, edges, edge_weights, node_weights):
        self.edges = edges
        self.edge_weights = edge_weights
        self.node_weights = node_weights
        self.n_parts = 2
        self.margin = 5.0
        self.coarsening_ratio = 3.0
        self.solutions = []
        self.solution_values = []
        self.initial_solutions = []

    def check(self):
        assert len(self.solutions) == len(self.solution_values)
        for s in self.solutions:
            assert len(s) == len(self.node_weights)
        for s in self.initial_solutions:
            assert len(s) == len(self.node_weights)

    def run(self):
        print ("Partitioner: " + str(len(self.node_weights)) + " nodes , " + str(len(self.edge_weights)) + " edges")
        target = len(self.node_weights) / self.coarsening_ratio
        max_run = int(2 * math.log(target, self.n_parts))
        self.n_merged = 1
        run = 0
        for i in range(max_run):
            run += 1
            self.run_once()
            self.infer_coarsening()
            if self.n_merged >= target:
                break
        self.sort_solutions()
        print ("Partitioner: " + str(run) + " runs")
        print ("Partitioner: coarsening ratio: " + str(len(self.node_weights) / self.n_merged) )
        print ("Partitioner: average cost: " + str(sum(self.solution_values) / len(self.solution_values)) )
        print ("Solutions: ", self.solution_values)
        sys.stdout.flush()
        if self.n_merged > 100:
            self.call_recursive()

    def call_recursive(self):
        next_level = self.build_recursive()
        next_level.check()
        next_level.run()

    def sort_solutions(self):
        sols = [a for a in zip(self.solution_values, self.solutions)]
        sols.sort()
        self.solution_values, self.solutions = zip(*sols)

    def build_recursive(self):
        # Build a coarsened graph
        coarsened_edges = []
        coarsened_edge_weights = []
        edge2weight = dict()
        for i in range(len(self.edges)):
            coarsened_edge = tuple(set([self.node2merged[p] for p in self.edges[i]]))
            edge2weight[coarsened_edge] = edge2weight.get(coarsened_edge, 0) + self.edge_weights[i]
        for edge, weight in edge2weight.items():
            coarsened_edges.append(edge)
            coarsened_edge_weights.append(weight)
        # Build the coarsened node weights
        coarsened_node_weights = [0 for i in range(self.n_merged)]
        for i in range(len(self.node_weights)):
            coarsened_node_weights[self.node2merged[i]] += self.node_weights[i]
        # Build the coarsened solutions
        coarsened_solutions = []
        for sol in self.solutions:
            coarsened = [0 for i in range(self.n_merged)]
            for i in range(len(self.node_weights)):
                coarsened[self.node2merged[i]] = sol[i]
            coarsened_solutions.append(coarsened)
        next_level = MultilevelSolver(coarsened_edges, coarsened_edge_weights, coarsened_node_weights)
        next_level.initial_solutions = coarsened_solutions
        return next_level

    def run_once(self):
        seed = len(self.solutions)
        builder = ModelBuilder(self.edges, self.edge_weights, self.node_weights)
        if len(self.initial_solutions) > 0:
            builder.initial_solution = self.initial_solutions.pop(0)
        builder.check()
        builder.build()
        builder.solve(seed=seed)
        self.solutions.append(builder.solution)
        self.solution_values.append(builder.objective_value)

    def infer_coarsening(self):
        node_placements = []
        for i in range(len(self.node_weights)):
            node_placements.append(tuple([s[i] for s in self.solutions]))
        self.node2merged = []
        place2merged = dict()
        merged_index = 0
        for i in range(len(self.node_weights)):
            pl = node_placements[i]
            if pl in place2merged:
                self.node2merged.append(place2merged[pl])
            else:
                place2merged[pl] = merged_index
                self.node2merged.append(merged_index)
                merged_index += 1
        self.n_merged = merged_index

if len(sys.argv) < 2:
    print("Usage: ls_partitioning.py graph.hgr")
    sys.exit(1)

graph_file_name = sys.argv[1]
random.seed(1)

edges, edge_weights, node_weights = read_graph(graph_file_name)
check_graph(edges, edge_weights, node_weights)

solver = MultilevelSolver (edges, edge_weights, node_weights)
solver.run()

#for i in range(10):
#    builder = ModelBuilder(edges, edge_weights, node_weights)
#    builder.build()
#    builder.solve(seed=i)

