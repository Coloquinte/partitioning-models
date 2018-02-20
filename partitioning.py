#!/usr/bin/python3

import localsolver
import sys
import random
import math

class Options():
    def __init__(self):
        # Partitions
        self.nb_parts = 2
        self.margin = 5.0
        # Cost function
        self.replication = False
        self.cost = "cut"
        # Multilevel
        self.coarsening_ratio = 3.0
        self.nb_cycles = 10
        # Solver
        self.verbosity = 0
        self.time_limit = None
        self.iteration_limit = None

class Graph(object):
    """
    Representation of the hypergraph
    """
    def __init__(self):
        self.edges = []
        self.edge_weights = []
        self.node_weights = []
        
    @staticmethod
    def read(graph_file_name):
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
        ret = Graph()
        ret.edges = edges
        ret.edge_weights = edge_weights
        ret.node_weights = node_weights
        return ret

    def nb_nodes (self):
        return len(self.node_weights)

    def nb_edges (self):
        return len(self.edge_weights)

    def check(self):
        for e in self.edges:
            for p in e:
                assert p >= 0 and p < self.nb_nodes()
        for w in self.edge_weights:
            assert w >= 0
        for w in self.node_weights:
            assert w >= 0
        assert len(self.edge_weights) == len(self.edges)

class ModelBuilder(object):
    """
    Create and solve a LocalSolver model
    """
    def __init__(self, graph, options):
        self.graph = graph 
        self.ls = localsolver.LocalSolver()
        self.model = self.ls.get_model()
        self.options = options
        # Constraints + starting point
        self.coarsening = []
        # Decisions
        self.node_placement = None
        # Debug
        self.edge_costs = []

    def check(self):
        self.graph.check()
        if self.node_placement != None:
            assert len(self.node_placement) == self.graph.nb_nodes()
    
    def check_solution(self, solution):
        if solution == None:
            return
        assert len(solution) == self.graph.nb_nodes()
        for s in solution:
            assert type(s) is tuple
            assert len(s) == self.options.nb_parts

    def solve(self, seed=0):
        self.ls.get_param().set_seed(seed)
        self.ls.solve()
        self.check()

    def build(self):
        self.check()
        self.build_node_variables()
        self.build_node_constraints()
        self.build_cost()
        self.apply_coarsening()
        # Close and apply parameters
        self.ls.get_param().set_verbosity(self.options.verbosity)
        self.model.close()
        self.ls.get_param().set_nb_threads(1)
        phase = self.ls.create_phase()
        if self.options.iteration_limit != None:
            phase.iteration_limit = self.options.iteration_limit
        if self.options.time_limit != None:
            phase.time_limit = self.options.time_limit
        self.check()

    def solution(self):
        s = []
        for placement in self.node_placement:
            s.append(tuple([b.value for b in placement]))
        return s

    def objective_value(self):
        return self.model.get_objective(0).value

    def build_cost(self):
        if self.options.replication:
            if self.options.cost == "cut":
                self.build_cut_cost_with_replication()
                #self.build_cut_cost()
            else:
                self.build_degree_cost_with_replication()
        else:
            if self.options.cost == "cut":
                self.build_cut_cost()
            else:
                self.build_degree_cost()

    def build_node_variables(self):
        self.node_placement = []
        for w in range(self.graph.nb_nodes()):
            places = [self.model.bool() for i in range(self.options.nb_parts)]
            self.node_placement.append(places)
            if self.options.replication:
                self.model.add_constraint(self.model.sum(places) >= 1)
            else:
                self.model.add_constraint(self.model.sum(places) == 1)

    def init_placement(self, solution=None):
        if solution == None:
            solution = []
            for i in range(self.graph.nb_nodes()):
                vals = [ False for p in range(self.options.nb_parts)]
                pos = random.randrange(self.options.nb_parts)
                vals[pos] = True
                solution.append(tuple(vals))
        self.check_solution(solution)
        for variables, values in zip(self.node_placement, solution):
            for variable, value in zip(variables, values):
                variable.set_value(value)

    def build_node_constraints(self):
        tot_weight = sum(self.graph.node_weights)
        self.weight_per_part = int(tot_weight * (1.0 + self.options.margin / 100.0) / self.options.nb_parts)
        for part in range(self.options.nb_parts):
            self.build_node_constraint(part)

    def build_node_constraint(self, part):
        weights_on_part = []
        for i in range(self.graph.nb_nodes()):
            decision = self.graph.node_weights[i]
            weight = self.node_placement[i][part]
            weights_on_part.append(decision * weight)
        self.model.add_constraint(self.model.sum(weights_on_part) < self.weight_per_part)

    def build_edge_degree(self, e):
        pins = self.graph.edges[e]
        counters = []
        for j in range(self.options.nb_parts):
            counters.append(self.model.sum([self.node_placement[p][j] for p in pins]))
        occupied = [c != 0 for c in counters]
        return self.model.sum(occupied)

    def build_edge_degree_without_source(self, e):
        """Cost function helper for node replication"""
        pins = self.graph.edges[e]
        if len(pins) <= 1:
            return 0
        sink_pins = pins[1:]
        source_pin = pins[0]
        counters = []
        for j in range(self.options.nb_parts):
            counters.append(self.model.sum([self.node_placement[p][j] for p in sink_pins]))
        source_present = []
        for j in range(self.options.nb_parts):
            source_present.append(self.node_placement[source_pin][j])
        m = self.model
        no_source = [m.and_(c !=  0, m.not_(s)) for c, s in zip(counters, source_present)]
        return self.model.sum(no_source)

    def apply_coarsening(self):
        for merged_group in self.coarsening:
            if len(merged_group) <= 1:
                continue
            n1 = merged_group[0]
            for n2 in merged_group[1:]:
                for a, b in zip(self.node_placement[n1], self.node_placement[n2]):
                    self.model.add_constraint(a == b)

    def build_cut_cost(self):
        self.edge_costs = [self.graph.edge_weights[i] * (self.build_edge_degree(i) >= 2) for i in range(self.graph.nb_edges())]
        cut = self.model.sum(self.costs)
        self.model.minimize(cut)

    def build_degree_cost(self):
        self.edge_costs = [self.graph.edge_weights[i] * self.build_edge_degree(i) for i in range(self.graph.nb_edges())]
        sum_degrees = self.model.sum(self.edge_costs) - sum(self.graph.edge_weights)
        self.model.minimize(sum_degrees)

    def build_cut_cost_with_replication(self):
        self.edge_costs = [self.graph.edge_weights[i] * (self.build_edge_degree_without_source(i) >= 1) for i in range(self.graph.nb_edges())]
        cut = self.model.sum(self.edge_costs)
        self.model.minimize(cut)

    def build_degree_cost_with_replication(self):
        self.edge_costs = [self.graph.edge_weights[i] * self.build_edge_degree_without_source(i) for i in range(self.graph.nb_edges())]
        cut = self.model.sum(self.edge_costs)
        self.model.minimize(cut)


class MultilevelSolver(object):
    """
    Use a search-driven coarsening algorithm to solve a partitioning problem
    """
    def __init__(self, graph, options):
        self.graph = graph
        self.options = options
        self.solutions = []
        self.coarsenings = []
        self.builder = None

    def solve(self):
        for i in range(options.nb_cycles):
            print ("Starting cycle #" + str(i+1))
            self.solve_recursive()

    def solve_recursive(self):
        while self.current_nodes() > 100:
            self.optimize_up()
            print ("Coarsening")
        while len(self.coarsenings) > 0:
            self.optimize_down()
            print ("Uncoarsening")
            self.coarsenings.pop()

    def optimize_up(self):
        solutions = self.solutions
        self.solutions = []
        self.objective_values = []
        coarsening = None
        for i in range(self.max_run()):
            builder = self.get_builder()
            solution = None
            if len(solutions) > 0:
                solution = solutions.pop(0)
            builder.init_placement(solution)
            builder.solve(seed=i)
            self.print_point()
            self.solutions.append(builder.solution())
            self.objective_values.append(builder.objective_value())
            coarsening = self.compute_coarsening()
            if len(coarsening) >= self.target_nodes():
                break
        print()
        self.sort_solutions()
        self.display()
        self.coarsenings.append(coarsening)

    def optimize_down(self):
        solutions = self.solutions
        self.solutions = []
        self.objective_values = []
        for i, solution in enumerate(solutions):
            builder = self.get_builder()
            builder.init_placement(solution)
            builder.solve(i)
            self.print_point()
            self.solutions.append(builder.solution())
            self.objective_values.append(builder.objective_value())
        print()
        self.sort_solutions()
        self.display()

    def print_point(self):
        print(".", end='')
        sys.stdout.flush()

    def display(self):
        print ("Clusters: " + str(self.current_nodes()) + ", nodes: " + str(self.graph.nb_nodes()) + ", edges: " + str(self.graph.nb_edges()))
        print (str(len(self.objective_values)) + " solutions: " + str(self.objective_values))
        print ("Average: " + str(sum(self.objective_values)/len(self.objective_values)))

    def get_builder(self):
        builder = ModelBuilder(self.graph, self.options)
        if len(self.coarsenings) > 0:
            builder.coarsening = self.coarsenings[-1]
        builder.build()
        return builder

    def compute_coarsening(self):
        placement2nodes = dict()
        for node in range(self.graph.nb_nodes()):
            placement = tuple([solution[node] for solution in self.solutions])
            placement2nodes.setdefault(placement, []).append(node)
        coarsening =  list(placement2nodes.values())
        return coarsening

    def sort_solutions(self):
        sols = [a for a in zip(self.objective_values, self.solutions)]
        sols.sort()
        self.objective_values = []
        self.solutions = []
        for o, s in sols:
            self.objective_values.append(o)
            self.solutions.append(s)

    def current_nodes(self):
        if len(self.coarsenings) == 0:
            return self.graph.nb_nodes()
        else:
            return len(self.coarsenings[-1])

    def target_nodes(self):
        clusters = self.current_nodes()
        return int(clusters / self.options.coarsening_ratio)

    def max_run(self):
        return int(2 * math.log(self.target_nodes(), 2))

if len(sys.argv) < 2:
    print("Usage: ls_partitioning.py graph.hgr margin")
    sys.exit(1)

graph_file_name = sys.argv[1]
if len(sys.argv) >= 4:
    random.seed(int(sys.argv[3]))
else:
    random.seed(1)

graph = Graph.read(graph_file_name)
graph.check()

options = Options()
if len(sys.argv) >= 3:
    options.margin = float(sys.argv[2])
options.iteration_limit = 50 * graph.nb_nodes()
options.time_limit = 10
options.replication = True
#options.verbosity=1

solver = MultilevelSolver(graph, options)
solver.solve()
#builder = ModelBuilder(graph, options)
#builder.build()
#builder.solve(seed)


