import heapq

def new_object(obj):
    return obj.__new__(obj)

__version__ = 1.0

cdef class MassHeap(object):
    cdef public  list contents
    cdef public MassCarrier carrier

    def __init__(self, items):
        self.carrier = MassCarrier(0.0)
        self.contents = []
        self.contents.extend(MassWrapper(i) for i in items)
        heapq.heapify(self.contents)
        self._sort()

    def __reduce__(self):
        return type(self), (list(self),), self.__getstate__()


    def __getstate__(self):
        d = {}
        d['items'] = list(self)
        return d

    def __setstate__(self, d):
        self.contents = []
        self.carrier = MassCarrier(0.0)
        self.contents.extend(MassWrapper(i) for i in d['items'])

    def __iter__(self):
        for x in self.contents:
            yield x.obj

    def _sort(self):
        self.contents.sort(key=lambda x: x.mass)

    def get_lower_than(self, float mass):
        cdef MassCarrier dummy = self.carrier
        dummy.mass = mass
        cdef MassWrapper item

        for item in self.contents:
            if item < dummy:
                yield item.obj
            else:
                break

    def __getitem__(self, index):
        return self.contents[index]

cdef class MassCarrier:
    cdef public  float mass
    def __init__(self, mass):
        self.mass = mass

    def __repr__(self):
        return "<Mass %f>" % self.mass

cdef class MassWrapper(MassCarrier):
    cdef public object obj

    def __init__(self, obj):
         self.obj = obj
         self.mass = obj.mass

    def __richcmp__(MassCarrier self, MassCarrier other, int opt):
        cdef float mmass = self.mass
        cdef float omass = other.mass
        if opt == 0:
            return mmass < omass
        elif opt == 2:
            return mmass == omass
        elif opt == 4:
            return mmass > omass
        elif opt == 1:
            return mmass <= omass
        elif opt == 5:
            return mmass >= omass
        elif opt == 3:
            return mmass != omass
