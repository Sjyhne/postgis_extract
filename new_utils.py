from data import Point, Edge
import networkx as nx
from shapely.geometry import Polygon, LineString
from shapely.validation import explain_validity

import numpy as np
from scipy.spatial import KDTree

def merge_close_nodes(edges, tolerance=0.01):
    # Extract all unique nodes
    unique_nodes = list({tuple(node) for edge in edges for node in [edge.start, edge.end]})
    kdtree = KDTree(unique_nodes)
    
    merged_nodes = {}
    for node in unique_nodes:
        idx = kdtree.query_ball_point(node, tolerance)
        merged_node = tuple(np.mean([unique_nodes[i] for i in idx], axis=0))
        for i in idx:
            merged_nodes[unique_nodes[i]] = merged_node
    
    merged_edges = []
    for edge in edges:
        merged_start = merged_nodes[tuple(edge.start)]
        merged_end = merged_nodes[tuple(edge.end)]
        merged_edges.append(Edge(Point(*merged_start), Point(*merged_end)))
    
    return merged_edges

def extract_unique_edges(data):
    unique_edges = set()
    for feature in data:
        if feature['type'] == 'MultiLineString':
            for linestring in feature['coordinates']:
                for i in range(len(linestring) - 1):
                    edge = Edge(linestring[i], linestring[i + 1])
                    unique_edges.add(edge)
    return unique_edges

def merge_collinear_edges(edges, tolerance_distance=0.001):
    merged_edges = set()
    edges = list(edges)
    while edges:
        edge = edges.pop()
        merged = False
        for other_edge in edges:
            if edge.is_collinear(other_edge):
                # Check distance between the ends of the edges
                if (LineString([edge.end, other_edge.start]).length <= tolerance_distance or
                    LineString([edge.start, other_edge.end]).length <= tolerance_distance):
                    new_edge = edge.merge_with(other_edge)
                    edges.remove(other_edge)
                    edges.append(new_edge)
                    merged = True
                    break
        if not merged:
            merged_edges.add(edge)
    return merged_edges

def build_graph(edges):
    graph = nx.Graph()  # Use an undirected graph
    for edge in edges:
        graph.add_edge(edge.start, edge.end)
    return graph


def find_unique_polygons(graph):
    cycles = list(nx.cycle_basis(graph))  # Ensure graph is undirected
    unique_polygons = []

    for cycle in cycles:
        if len(cycle) >= 3:
            polygon = [Point(node[0], node[1], node[2]) for node in cycle]
            shapely_polygon = Polygon(polygon)
            if shapely_polygon.is_valid:
                unique_polygons.append(shapely_polygon)

    return unique_polygons

def convert_to_geojson_feature(polygon):
    return {
        "type": "Feature",
        "geometry": polygon.__geo_interface__,
        "properties": {}
    }

def validate_polygons(polygons):
    valid_polygons = []
    for polygon in polygons:
        if len(polygon.exterior.coords) >= 4 and polygon.is_valid:
            valid_polygons.append(polygon)
    return valid_polygons

def subtract_smaller_polygons_from_larger(polygons):
    modified_polygons = []
    for i, poly in enumerate(polygons):
        modified_poly = poly
        for j, other_poly in enumerate(polygons):
            if i != j and poly.contains(other_poly):
                modified_poly = modified_poly.difference(other_poly)
        modified_polygons.append(modified_poly)
    return modified_polygons

def filter_redundant_polygons(polygons):
    unique_polygons = []
    for i, poly in enumerate(polygons):
        is_composed = False
        for j, other_poly in enumerate(polygons):
            if i != j and other_poly.contains(poly):
                difference = other_poly.difference(poly)
                if difference.is_empty:
                    is_composed = True
                    break
        if not is_composed:
            unique_polygons.append(poly)
    return unique_polygons


def filter_duplicate_polygons(polygons):
    unique_polygons = []
    
    for poly in polygons:
        if poly not in unique_polygons:
            unique_polygons.append(poly)
            
    return unique_polygons