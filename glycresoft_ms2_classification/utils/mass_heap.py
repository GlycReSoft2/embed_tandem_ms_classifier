import heapq
try:
    from .cmass_heap import MassHeap as CMassHeap
    use_cython = True
except:
    use_cython = False


class PMassHeap(object):
    '''
        Defines a heap based on the mass of each element in its contents.

        If its elements do not possess a `mass` attribute, then a key function
        must be provided to calculate the mass for location
    '''
    def __init__(self, items=None, key=None):
        self.contents = []
        self.carrier = MassCarrier(0.0)
        if items is not None:
            self.contents.extend(MassWrapper(i, key) for i in items)
            heapq.heapify(self.contents)
            self._sort()

    def add(self, item, key=None):
        item = MassWrapper(item, key)
        heapq.heappush(self.contents, item)

    def extend(self, items, key=None):
        for item in items:
            self.add(item, key)

    def __iter__(self):
        gen = iter(self.contents)
        for i in gen:
            yield i.object

    def descending(self):
        gen = reversed(self.contents)
        for i in gen:
            yield i.object

    ascending = __iter__

    def _sort(self):
        self.contents.sort()

    def get_lower_than(self, mass):
        dummy = self.carrier
        dummy.mass = mass
        for item in self.contents:
            if item < dummy:
                yield item.object
            else:
                break

    def __repr__(self):
        return repr(self.contents)

    def __getitem__(self, index):
        return self.contents[index]


class MassCarrier(object):
    def __init__(self, mass):
        self.mass = mass


class MassWrapper(object):
    def __init__(self, object, key=None):
        self.object = object
        self.mass = object.mass if key is None else key(object)

    def __lt__(self, other):
        return self.mass < other.mass

    def __gt__(self, other):
        return self.mass > other.mass

    def __repr__(self):
        return "MassWrapper({0})[{1}]".format(str(self.object), self.mass)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        '''Comparison should be based on comparing
        two wrappers with object attributes'''
        try:
            return self.object == other.object
        except:
            return False

    def __ne__(self, other):
        try:
            return self.object != other.object
        except:
            return True

if use_cython:
    MassHeap = CMassHeap
else:
    MassHeap = PMassHeap
