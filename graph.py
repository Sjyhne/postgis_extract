import math

from shapely.geometry import Point, LineString
import rtree
from shapely import Polygon
import uuid

import numpy as np

def point_to_line_distance(point, line_start, line_end):
    # Convert points to numpy arrays for vector operations
    point = np.array(point)
    line_start = np.array(line_start)
    line_end = np.array(line_end)
    
    # Calculate the vector from the line start to the point and from the line start to the line end
    line_vec = line_end - line_start
    point_vec = point - line_start
    
    # Calculate the projection of the point vector onto the line vector
    line_len_squared = np.dot(line_vec, line_vec)
    point_projection = np.dot(point_vec, line_vec) / line_len_squared
    
    # Ensure the projection point is on the line segment
    if point_projection < 0.0:
        closest_point = line_start
    elif point_projection > 1.0:
        closest_point = line_end
    else:
        closest_point = line_start + point_projection * line_vec
    
    # Return the distance from the point to the closest point on the line segment
    return np.linalg.norm(point - closest_point)

# Create an R-tree index
index = rtree.index.Index()

def distance(a, b):
    """Calculate the Euclidean distance between two points."""
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5

def is_collinear(node1, node2, node3, epsilon=1e-6):
    x1, y1 = node1.x, node1.y
    x2, y2 = node2.x, node2.y
    x3, y3 = node3.x, node3.y
    determinant = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
    return abs(determinant) < epsilon

def point_line_distance(point, start, end):
    """Calculate the perpendicular distance from a point to a line segment."""
    if (start[0] == end[0]) and (start[1] == end[1]):
        return math.sqrt((point[0] - start[0])**2 + (point[1] - start[1])**2)

    # Calculate the line equation coefficients
    A = end[1] - start[1]
    B = start[0] - end[0]
    C = A*start[0] + B*start[1]

    # Calculate the distance
    distance = abs(A*point[0] + B*point[1] - C) / math.sqrt(A**2 + B**2)
    return distance

def cross_product(a, b, c):
    """Calculate the cross product AB x AC."""
    ab = (b[0] - a[0], b[1] - a[1])
    ac = (c[0] - a[0], c[1] - a[1])
    return ab[0] * ac[1] - ab[1] * ac[0]

def is_almost_on_line(a, b, c, tolerance=0.3):
    """Check if point c is almost on the line segment between a and b."""
    cp = cross_product(a, b, c)
    if abs(cp) > tolerance:
        return False  # c is not on the line if the cross product is not close to zero
    
    dot_product = (c[0] - a[0]) * (b[0] - a[0]) + (c[1] - a[1]) * (b[1] - a[1])
    squared_length_ab = (b[0] - a[0])**2 + (b[1] - a[1])**2
    return 0 <= dot_product <= squared_length_ab 


def is_point_between(a, b, c, margin=1.0):
    """
    Check if point c is near the line segment ab within a specified margin.
    """
    # Unpack points
    xa, ya = a
    xb, yb = b
    xc, yc = c
    
    # Calculate the numerator of the distance formula
    numerator = abs((yb - ya) * xc - (xb - xa) * yc + xb * ya - yb * xa)
    # Calculate the denominator of the distance formula
    denominator = ((yb - ya)**2 + (xb - xa)**2)**0.5
    # Calculate the perpendicular distance
    distance = numerator / denominator if denominator else float('inf')
    
    # Check if distance is within the margin
    if distance <= margin:
        # Ensure C projects onto the segment AB, not just the line
        dot_product = (xc - xa) * (xb - xa) + (yc - ya) * (yb - ya)
        len_sq_ab = (xb - xa)**2 + (yb - ya)**2
        # Check if C's projection lies between A and B
        return 0 <= dot_product <= len_sq_ab
    return False

def is_point_on_line_segment(p1, p2, p, epsilon=1e-6):
    # Cross product to check if the points are collinear
    cross_product = (p[1] - p1[1]) * (p2[0] - p1[0]) - (p[0] - p1[0]) * (p2[1] - p1[1])
    if abs(cross_product) > epsilon:
        return False  # Not collinear

    # Dot product and squared length to check if p is between p1 and p2
    dot_product = (p[0] - p1[0]) * (p2[0] - p1[0]) + (p[1] - p1[1]) * (p2[1] - p1[1])
    squared_length_p1_p2 = (p2[0] - p1[0])**2 + (p2[1] - p1[1])**2
    return 0 <= dot_product <= squared_length_p1_p2


class Node:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z if z >= 0 else 0
        self.idx = str(uuid.uuid4())
        self.neighbors = set()  # Using a set to avoid duplicate neighbors

    def add_neighbor(self, neighbor):
        self.neighbors.add(neighbor)

    def __eq__(self, other):
        return (round(self.x, 3), round(self.y, 3)) == (round(other.x, 3), round(other.y, 3))

    def __hash__(self):
        return hash((self.x, self.y))

    def __repr__(self):
        return f"Node({self.x}, {self.y}, {self.z}, Neighbors={len(self.neighbors)})"

class Graph:
    def __init__(self):
        self.nodes = set()

    def add_or_get_node(self, tmp_node):
        for node in self.nodes:
            if node == tmp_node:
                return node
        self.nodes.add(tmp_node)
        return tmp_node

    def connect_nodes(self, node, other_node):
        node1 = self.add_or_get_node(node)
        node2 = self.add_or_get_node(other_node)
        node1.add_neighbor(node2)
        node2.add_neighbor(node1)

    def point_to_line_distance(self, point, line_start, line_end):
        point = np.array((point.x, point.y))
        line_start = np.array((line_start.x, line_start.y))
        line_end = np.array((line_end.x, line_end.y))
        line_vec = line_end - line_start
        point_vec = point - line_start
        line_len_squared = np.dot(line_vec, line_vec)
        point_projection = np.dot(point_vec, line_vec) / line_len_squared
        if point_projection < 0.0 or point_projection > 1.0:
            return np.inf  # The point is not within the line segment
        closest_point = line_start + point_projection * line_vec
        return np.linalg.norm(point - closest_point)

    def adjust_nodes_close_to_edges(self, threshold):
        for node in list(self.nodes):  # Convert to list to avoid modification during iteration
            for neighbor in list(node.neighbors):
                for potential_close_node in self.nodes - {node, neighbor}:
                    if self.point_to_line_distance(potential_close_node, node, neighbor) < threshold:
                        # Connect the potential_close_node with node and neighbor, making it part of the edge
                        self.connect_nodes(potential_close_node, node)
                        self.connect_nodes(potential_close_node, neighbor)
                        # Optionally, remove the direct connection between node and neighbor if necessary
                        # This depends on whether you want to "split" the edge or just add the node as an intermediate point
                        node.neighbors.discard(neighbor)
                        neighbor.neighbors.discard(node)
    
    def remove_unnecessary_nodes(self, collinearity_threshold=0.1):
        removed_nodes = set()
        for node in list(self.nodes):  # Convert to list to avoid modification during iteration
            if node in removed_nodes:
                continue
            if len(node.neighbors) == 2:  # Node has exactly two neighbors
                neighbors = list(node.neighbors)
                n1, n2 = neighbors[0], neighbors[1]
                if self.is_almost_collinear(n1, node, n2, collinearity_threshold):
                    # Connect the two neighbors directly, bypassing the current node
                    self.connect_nodes(n1, n2)
                    # Remove the current node from its neighbors
                    n1.neighbors.discard(node)
                    n2.neighbors.discard(node)
                    # Finally, remove the node from the graph
                    self.nodes.remove(node)
                    removed_nodes.add(node)

    def remove_almost_collinear_nodes(self):
        nodes_to_remove = []
        for node in self.nodes:
            if len(node.neighbors) == 2:  # Node has exactly two neighbors
                n1, n2 = list(node.neighbors)
                if is_collinear(n1, node, n2):
                    nodes_to_remove.append(node)
        
        # Remove identified nodes and adjust connections
        for node in nodes_to_remove:
            n1, n2 = list(node.neighbors)
            # Connect the two neighbors to each other if they aren't already
            n1.neighbors.add(n2)
            n2.neighbors.add(n1)
            # Remove the current node from its neighbors
            n1.neighbors.remove(node)
            n2.neighbors.remove(node)
            # Remove the node from the graph
            self.nodes.remove(node)

    def is_almost_collinear(self, n1, n2, n3, threshold):
        # Reuse the point_to_line_distance function to check if n2 is almost on the line between n1 and n3
        distance = self.point_to_line_distance(n2, n1, n3)
        return distance < threshold

    def find_smallest_polygons(self):
        smallest_polygons = []
        already_searched = set()

        for start_node in self.nodes:
            if start_node in already_searched:
                continue
            
            result = {'nodes': [], 'area': float('inf')}
            visited = set()
            path = []
            self._find_smallest_polygon_dfs(start_node, start_node, path, visited, result, first_call=True, already_searched=already_searched)
            if result['nodes']:
                smallest_polygons.append(result)
            already_searched.add(start_node)
        
        return smallest_polygons

    def _find_smallest_polygon_dfs(self, current_node, start_node, path, visited, result, first_call=False, already_searched=set()):
        MAX_PATH_LENGTH = 30
        if len(path) > MAX_PATH_LENGTH:
            return

        if current_node in visited and not first_call:
            return

        path.append(current_node)
        if not first_call:
            visited.add(current_node)

        if len(path) > 3 and current_node == start_node and len(set(path[:-1])) == len(path) - 1:
            area = self._calculate_polygon_area(path)
            if area < result['area'] and area > 0:
                result['area'] = area
                result['nodes'] = path[:-1]

        else:
            for neighbor in current_node.neighbors - visited:
                self._find_smallest_polygon_dfs(neighbor, start_node, path.copy(), visited.copy(), result, False, already_searched)
        
        path.pop()
        if not first_call:
            visited.remove(current_node)

    def get_node_neighbors(self):
        node_neighbors = {}
        for node in self.nodes:
            neighbors = [neighbor for neighbor in node.neighbors]
            node_neighbors[f"{node.x},{node.y},{node.z}"] = [(neighbor.x, neighbor.y, neighbor.z) for neighbor in neighbors]
        return node_neighbors

    def _calculate_polygon_area(self, nodes):
        poly = Polygon([(node.x, node.y, node.z) for node in nodes])
        return poly.area

    def __repr__(self):
        return f"Graph({len(self.nodes)} nodes)"
    
    def validate_polygon(self, polygon_nodes):
        """
        Validates a polygon defined by a list of nodes using Shapely.
        :param polygon_nodes: List of nodes defining the polygon.
        :return: True if the polygon is valid, False otherwise.
        """
        # Convert nodes to a list of tuples (x, y)
        coords = [(node.x, node.y, node.z) for node in polygon_nodes]
        # Create a Shapely Polygon
        poly = Polygon(coords)
        # Return the validity
        return poly.is_valid
    
    def buffer_polygon(self, polygon_nodes):
        pass

    def to_geojson(self, polygons):
        """
        Convert a list of polygons into GeoJSON, validating each polygon.
        :param polygons: A list of polygons, where each polygon is represented by a list of nodes.
        """
        valid_multi_polygon_coords = []

        for polygon in polygons:
            # Validate the polygon
            if self.validate_polygon(polygon['nodes']):
                # Prepare coordinates for each polygon, ensuring closure
                coordinates = [[[node.x, node.y, node.z] for node in polygon['nodes'] + [polygon['nodes'][0]]]]
                valid_multi_polygon_coords.append(coordinates)

        # Construct the MultiPolygon feature only with valid polygons
        feature = {
            "type": "Feature",
            "properties": {}, # "building_box": building_box, "building_id": building_id},
            "geometry": {
                "type": "MultiPolygon",
                "coordinates": valid_multi_polygon_coords
            }
        }

        return feature