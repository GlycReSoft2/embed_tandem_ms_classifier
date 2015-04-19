import logging
from collections import defaultdict
from math import fabs

from components import masser, generate_component_set, unordered_combinations, SequenceCollection, SimpleFragment
from glycresoft_ms2_classification.structure import modification
from glycresoft_ms2_classification.structure.composition import Composition
from glycresoft_ms2_classification.utils.mass_heap import MassHeap

logging.basicConfig(level='DEBUG', format='%(asctime)s %(message)s', datefmt='%I:%M:%S %p')
logger = logging.getLogger()


precursor_mass_shift = Composition('H2O').mass
y_mass_shift = Composition('H2O').mass + Composition('H').mass - Composition('e').mass
b_mass_shift = Composition('H').mass - Composition('e').mass


def ppm_error(query, match):
    return (query - match)/match


def strings(graph):
    for seq in graph.all_sequences():
        for concrete in seq[1][1].to_sequence():
            yield concrete


def match(path, sequence):
    sequence = SimpleFragment(str(sequence)).to_sequence_positions()
    i = 0
    j = 0
    matched = 0
    current_chunk = set(path[j].link_sequence)
    while(i < len(sequence)):
        if sequence[i] in current_chunk:
            current_chunk.discard(sequence[i])
            matched += 1
            if len(current_chunk) == 0:
                j += 1
                current_chunk = set(path[j].link_sequence)
    return matched


class Node(object):
    def __init__(self, mass, composition=None, kind=None, graph=None):
        self.mass = mass
        self.kind = kind or []
        self.composition = composition or[]
        self.graph = graph

    def __hash__(self):
        return hash((self.mass))

    def __eq__(self, other):
        return (self.mass == other.mass) and (self.kind == other.kind)

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return "<Node {} {}>".format(self.mass, list(zip(self.composition, self.kind)))


class Edge(object):
    def __init__(self, from_terminus, to_terminus, base_mass, link_sequence, kind=''):
        self.from_terminus = from_terminus
        self.to_terminus = to_terminus
        self.base_mass = base_mass
        self.link_sequence = link_sequence
        self.kind = kind

    def __hash__(self):
        return hash((self.from_terminus, self.to_terminus, self.link_sequence))

    def __eq__(self, other):
        return self.from_terminus == other.from_terminus and\
            self.to_terminus == other.to_terminus and\
            self.link_sequence == other.link_sequence and\
            self.kind == other.kind

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        rep = "<{from_}-{link}-{to} {kind}>".format(
            from_=self.from_terminus, link=self.link_sequence,
            to=self.to_terminus, kind=self.kind)
        return rep

    def __len__(self):
        return len(self.link_sequence)


class Graph(object):
    def __init__(self, nodes, constant_modifications=None, variable_modifications=None):
        self.nodes = MassHeap(Node(s.mass, graph=self) for s in nodes)
        self.edges = defaultdict(list)
        self.constant_modifications = constant_modifications or []
        self.variable_modifications = variable_modifications or []
        self.variable_modifications += [modification.Modification("HexNAc").rule]
        self.parts = MassHeap(generate_component_set(
            self.constant_modifications,
            self.variable_modifications))
        self.long_parts = {1: self.parts}
        self.node_map = {}

    def __iter__(self):
        return iter(self.nodes)

    def process(self, node, parts=None):
        if parts is None:
            parts = self.parts
        for extent in self.nodes.get_higher_than(node.mass):
            gap_mass = extent.mass - node.mass
            for part in parts.get_lower_than(gap_mass + precursor_mass_shift + 1):
                if fabs(ppm_error(gap_mass + y_mass_shift, part.mass + y_mass_shift)) <= 2e-5:
                    self.edges[frozenset((node, extent))].append(Edge(
                            node, extent, node.mass, link_sequence=part)
                    )

    def process_all(self, length=1):
        parts = self.long_parts.get(length, None)
        if parts is None:
            self.build_unordered_sequences(length)
            parts = self.long_parts[length]
        for node in self:
            self.process(node, parts)
        self.build_node_map()

    def build_node_map(self):
        self.node_map = defaultdict(list)
        for node in self:
            for pair in [pair for pair in self.edges if node in pair]:
                if node == min(pair, key=masser):
                    self.node_map[node].extend(self.edges[pair])

    def roots(self):
        roots = set(self.node_map)
        for node, edges in self.node_map.items():
            for edge in edges:
                roots.discard(edge.to_terminus)
            if len(edges) == 0:
                roots.discard(node)
        return list(roots)

    def build_unordered_sequences(self, n=2):
        self.long_parts[n] = MassHeap(list(unordered_combinations(self.parts, n)))

    def get_sequence(self, node):
        for path in self.traverse(node):
            yield (node, SequenceCollection(map(lambda x: x.link_sequence, path)), path[-1].to_terminus)

    def traverse(self, node):
        if len(self.node_map[node]) == 0:
            yield []
        else:
            for edge in self.node_map[node]:
                for path in self.traverse(edge.to_terminus):
                    yield [edge] + path

    def all_paths(self):
        for root in self.roots():
            for path in self.traverse(root):
                yield (sum(map(len, path)), path)

    def all_sequences(self):
        for root in self.roots():
            for seq in self.get_sequence(root):
                yield (len(seq[1]), seq)

    def identify_node(self, node, parts=None):
        if parts is None:
            parts = self.parts
        for part in parts:
            if fabs(ppm_error(node.mass, part.mass + y_mass_shift)) <= 2e-5:
                node.composition.append(part)
                node.kind.append('y')
            elif fabs(ppm_error(node.mass, part.mass + b_mass_shift)) <= 2e-5:
                node.composition.append(part)
                node.kind.append('b')
        return zip(node.composition, node.kind)
