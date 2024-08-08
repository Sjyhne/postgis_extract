from collections import defaultdict
from shapely.geometry import LineString

class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    
    def __iter__(self):
        return iter((self.x, self.y, self.z))
    
    def __eq__(self, other):
        return (self.x, self.y, self.z) == (other.x, other.y, other.z)
    
    def __hash__(self):
        return hash((self.x, self.y, self.z))
    
    def __repr__(self):
        return f"Point({self.x}, {self.y}, {self.z})"

class Edge:
    def __init__(self, start, end):
        self.start = tuple(start)
        self.end = tuple(end)
        self.normalize()

    def normalize(self):
        if self.start > self.end:
            self.start, self.end = self.end, self.start

    def is_collinear(self, other):
        # Check if two edges are collinear
        line1 = LineString([self.start, self.end])
        line2 = LineString([other.start, other.end])
        return line1.intersection(line2).equals(line1) or line1.intersection(line2).equals(line2)

    def merge_with(self, other):
        # Merge two collinear edges
        return Edge(self.start, other.end)

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def __hash__(self):
        return hash((self.start, self.end))

    def __repr__(self):
        return f"Edge({self.start}, {self.end})"


class Graph:
    def __init__(self):
        self.adjacency_list = defaultdict(set)
    
    def add_edge(self, edge):
        self.adjacency_list[edge.start].add(edge.end)
        self.adjacency_list[edge.end].add(edge.start)
    
    def __repr__(self):
        return f"Graph({dict(self.adjacency_list)})"