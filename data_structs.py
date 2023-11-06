from shapely.geometry import Point, Polygon as ShapelyPolygon

class Node:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.neighbours = list()
        self.edges = []  # List of edges connected to this node

    def __eq__(self, other):
        if other is None:
            return False
        return (self.x, self.y, self.z) == (other.x, other.y, other.z)
    
    def __hash__(self):
        return hash((self.x, self.y, self.z))

    def __repr__(self) -> str:
        return f"Node({self.x}, {self.y}, {self.z})"

    def add_edge(self, edge):
        if edge not in self.edges:
            self.edges.append(edge)

class Edge:
    def __init__(self, node1, node2):
        self.node1 = node1
        self.node2 = node2
        self.used = False
    
    def __repr__(self) -> str:
        return f"Edge({self.node1}, {self.node2})"
    
class Polygon:
    def __init__(self, nodes) -> None:
        self.nodes = nodes
    
    def calculate_area(self):
        if len(self.nodes) < 3:
            return 0.0
        
        area = 0.0
        for i in range(len(self.nodes)):
            x1, y1, _ = self.nodes[i].x, self.nodes[i].y, self.nodes[i].z
            x2, y2, _ = self.nodes[(i + 1) % len(self.nodes)].x, self.nodes[(i + 1) % len(self.nodes)].y, self.nodes[(i + 1) % len(self.nodes)].z
            area += x1 * y2 - x2 * y1
        
        self.area = abs(area) / 2.0

class Graph:
    def __init__(self):
        self.nodes = []
        self.edges = []

    def add_node(self, node):
        self.nodes.append(node)
        return node

    def connect_nodes(self, node1, node2):
        edge = Edge(node1, node2)
        self.edges.append(edge)
        node1.add_edge(edge)
        node2.add_edge(edge)
    
    def populate_node_neighbours(self):
        """Populate each node's neighbours using the Edge objects."""
        for edge in self.edges:
            # Add the nodes to each other's neighbours if not already added
            if edge.node2 not in edge.node1.neighbours:
                edge.node1.neighbours.append(edge.node2)
            if edge.node1 not in edge.node2.neighbours:
                edge.node2.neighbours.append(edge.node1)
    
    def remove_duplicate_nodes(self):
        unique_nodes = {}
        for node in self.nodes:
            key = (node.x, node.y, node.z)
            unique_nodes[key] = node

        # Update nodes list to unique nodes
        self.nodes = list(unique_nodes.values())

        # Reconstruct edges to ensure no invalid edges remain
        new_edges = set()  # Using a set to ensure no duplicate edges
        for edge in self.edges:
            key1, key2 = (edge.node1.x, edge.node1.y, edge.node1.z), (edge.node2.x, edge.node2.y, edge.node2.z)
            if key1 != key2:  # Ensure the edge connects two distinct nodes
                # Sort the keys to ensure no reverse duplicate edges
                sorted_key1, sorted_key2 = tuple(sorted([key1, key2]))
                new_edges.add((sorted_key1, sorted_key2))

        # Update edges list
        self.edges = [Edge(unique_nodes[key1], unique_nodes[key2]) for key1, key2 in new_edges]

        # Update node edges
        for node in self.nodes:
            node.edges = [edge for edge in self.edges if edge.node1 == node or edge.node2 == node]
    
    def is_simple_polygon(self, polygon):
        """
        Check if the given polygon is simple (no other nodes lie within it).
        """
        # Create the polygon object using Shapely
        poly = ShapelyPolygon([(node.x, node.y) for node in polygon])

        # Check each node in the graph
        for node in self.nodes:
            point = Point(node.x, node.y)
            # If the node is not on the polygon's boundary and it's inside the polygon, then it's not simple
            if not poly.boundary.contains(point) and poly.contains(point):
                return False

        return True
    
    def find_polygons(self):
        """Find all polygons in the graph."""
        visited = set()
        stack = []
        polygons = []

        def dfs(node, parent):
            if node in visited:
                # We've found a cycle. Try to build a polygon from the stack
                # up to the first occurrence of the node.
                if node not in stack:
                    return
                idx = stack.index(node)
                potential_polygon = stack[idx:]
                # Verify if this is a simple polygon
                polygons.append(potential_polygon)
                return
            visited.add(node)
            stack.append(node)
            for neighbor in node.neighbours:
                if neighbor != parent:
                    dfs(neighbor, node)
            stack.pop()

        for node in self.nodes:
            if node is None:
                continue
            if node not in visited:
                dfs(node, None)

        return polygons