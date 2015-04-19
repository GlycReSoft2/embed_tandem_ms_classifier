from collections import defaultdict
from itertools import chain
from glycresoft_ms2_classification.structure import sequence
Sequence = sequence.Sequence


def simple_mass(res, mods):
    mass = res.mass
    for mod in mods:
        mass += mod.mass
    return mass


def freeze_pair(pair):
    return (pair[0], tuple(pair[1]))


def freeze_sequence(seq):
    collect = [(res, tuple(mods)) for res, mods in seq]
    return tuple(collect)


def extend_node(node, part, gap=False):
    mass = node.mass
    mass += part.mass
    next_node = Node(part, mass=mass, parents={node.part: node},
                     total_gaps=node.total_gaps + gap,
                     current_gaps=node.current_gaps + 1 if gap else 0)
    return next_node


def get_first_leaf(subroot):
    for node in subroot.children.values():
        if len(node.children) == 0:
            return node
        else:
            return get_first_leaf(node)


def get_first_root(leaf):
    for node in leaf.parents.values():
        if len(node.parents) == 0:
            return node
        else:
            return get_first_root(node)


class SequencePosition(tuple):
    def __new__(self, residue, modifications):
        return tuple.__new__(self, (residue, modifications))

    @property
    def mass(self):
        mass = self[0].mass
        for mod in self[1]:
            mass += mod.mass
        return mass


class Node(object):
    def __init__(self, part, mass=None, index=0, parents=None, children=None,
                 current_gaps=0, total_gaps=0, graph=None):
        self.part = part
        self.parents = parents or {}
        self.children = children or {}
        self.mass = mass
        self.index = index
        self.graph = graph
        for node in self.parents.values():
            node.add_child(self)
        self.current_gaps = current_gaps
        self.total_gaps = total_gaps

    def to_sequence(self):
        if len(self.children) == 0:
            yield (self.part,)
        for child in self.children.values():
            for child_seq in chain.from_iterable(child.to_sequence()):
                yield (self.part,) + child_seq

    def to_path(self):
        if len(self.children) == 0:
            yield (self,)
        for child in self.children.values():
            for child_seq in child.to_path():
                yield (self,) + child_seq

    def add_parent(self, parent_node):
        self.parents[parent_node.part] = (parent_node)

    def add_child(self, child_node):
        self.children[(child_node.part)] = child_node
        if child_node.graph is None:
            child_node.graph = self.graph
        child_node.add_parent(self)

    def __repr__(self):
        rep = "<{} {} {}/{}>".format(self.part, self.mass, self.current_gaps, self.total_gaps)
        # for child in self.children.values():
        #     rep += "\n\t"
        #     child_parts = repr(child).split('\n')
        #     header = child_parts[0]
        #     rep += header + '\n'
        #     rep += '\n\t'.join(child_parts[1:])
        return rep

    def reverse(self):
        parents = self.parents
        children = self.children
        self.parents = children
        self.children = parents


class SequenceGraph(object):

    def __init__(self, roots=None):
        self.roots = roots or {}
        if not isinstance(self.roots, dict):
            root_dict = {}
            for root in self.roots:
                root_dict[root.part] = root
            self.roots = root_dict
        self.nodes = defaultdict(dict)

    def get(self, index, seq_pos):
        return self.nodes[index].get(seq_pos)

    def traverse(self, seq):
        if not isinstance(seq, Sequence):
            seq = Sequence(seq)
        seq = freeze_sequence(seq)
        gen = iter(seq)
        start = self.roots[gen.next()]
        layer = []
        next_layer = [[[], start]]
        for child in gen:
            layer = next_layer
            next_layer = []
            for path, node in layer:
                next_layer.append([path + [node], node.children[child]])
        return next_layer

    @classmethod
    def from_sequence(cls, seq, shift=0):
        if not isinstance(seq, Sequence):
            seq = Sequence(seq)
        seq = freeze_sequence(seq)
        gen = iter(seq)
        start = gen.next()
        mass = simple_mass(*start) + shift
        root = Node(start, mass)
        prev = root
        for pos in gen:
            mass += simple_mass(*pos)
            pos = Node(pos, mass)
            prev.add_child(pos)
            prev = pos
        return cls(roots={root.part: root})

    def to_sequence(self):
        for root in self.roots.values():
            for seq in root.to_sequence():
                yield Sequence.from_iterable(seq)

    @property
    def leaves(self):
        for root in self.roots.values():
            stack = [root]
            while len(stack) > 0:
                node = stack.pop()
                if len(node.children) == 0:
                    yield node
                else:
                    stack.extend(node.children.values())
