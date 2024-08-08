
import logging
from collections import defaultdict
from shapely.geometry import Polygon


# Set up logging configuration
logging.basicConfig(level=logging.INFO, format='%(message)s')

def clean_graph(graph):
    clean_graph = defaultdict(set)
    
    for node, connections in graph.items():
        node_str = str(node).strip("() ")
        for conn in connections:
            conn_str = str(conn).strip("() ")
            clean_graph[node_str].add(conn_str)
            clean_graph[conn_str].add(node_str)
    
    return {k: list(v) for k, v in clean_graph.items()}

def construct_polygons(graph):
    def find_polygon(start_node, graph, visited_edges):
        visited = set()
        stack = [(start_node, None)]
        polygon = []
        
        while stack:
            node, parent = stack.pop()
            if node not in visited:
                visited.add(node)
                polygon.append(node)
                for neighbor in graph[node]:
                    edge = tuple(sorted([node, neighbor]))
                    if neighbor != parent and edge not in visited_edges:
                        stack.append((neighbor, node))
                        visited_edges.add(edge)
        
        # Ensure the polygon is closed
        if polygon and polygon[0] != polygon[-1]:
            polygon.append(polygon[0])
        
        return polygon

    polygons = []
    visited_edges = set()
    
    for node in graph:
        polygon = find_polygon(node, graph, visited_edges)
        if polygon and is_valid_polygon(polygon):
            polygons.append(polygon)
    
    # Remove duplicate polygons
    unique_polygons = []
    seen_polygons = set()
    for polygon in polygons:
        if len(polygon) >= 4:  # Ensure at least 3 unique nodes + 1 to close the loop
            poly_tuple = tuple(sorted(polygon))
            if poly_tuple not in seen_polygons:
                seen_polygons.add(poly_tuple)
                unique_polygons.append(polygon)
    
    return unique_polygons

def is_valid_polygon(polygon):
    if len(polygon) < 4:  # A valid polygon must have at least 3 nodes plus a closing node
        return False
    
    poly_coords = [(float(n.split(',')[0]), float(n.split(',')[1])) for n in polygon]
    shapely_polygon = Polygon(poly_coords)
    
    # Check if the polygon is simple (edges do not intersect)
    if not shapely_polygon.is_simple:
        return False
    
    return True

def convert_to_geojson_feature(polygons):
    
    multipolygon = []
    
    for polygon in polygons:
        coordinates = []
        for node in polygon:
            lon, lat, alt = map(float, node.split(','))
            coordinates.append([lon, lat, alt])
        multipolygon.append([coordinates])
    
    feature = {
        "type": "Feature",
        "properties": {},
        "geometry": {
            "type": "MultiPolygon",
            "coordinates": multipolygon
        }
    }
    
    return feature


def generate_query_boxes(xmin, ymin, xmax, ymax, box_width, box_height):
    smaller_boxes = []

    current_xmin = xmin
    while current_xmin < xmax:
        current_ymin = ymin
        while current_ymin < ymax:
            # Calculate the top-right coordinates of the smaller box
            current_xmax = min(current_xmin + box_width, xmax)
            current_ymax = min(current_ymin + box_height, ymax)

            # Add the smaller box to the list
            smaller_boxes.append((current_xmin, current_ymin, current_xmax, current_ymax))

            # Move to the next box vertically
            current_ymin += box_height

        # Move to the next box horizontally
        current_xmin += box_width

    return smaller_boxes