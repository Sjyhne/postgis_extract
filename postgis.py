from postgisutils import get_cursor

import re
import json

from structures.point import Point
from structures.linestringz import LineStringZ

from collections import defaultdict


def parse_linestringz(wkt):
    """
    Parses a LINESTRING Z in WKT format and returns coordinates as a list of tuples.
    
    Args:
    - wkt (str): The LINESTRING Z in WKT format.
    
    Returns:
    - List[Tuple[float, float, float]]: A list of tuples containing the X, Y, and Z coordinates.
    """
    
    if type(wkt) == tuple:
        wkt = wkt[0]
    
    # Extract the SRID
    srid_match = re.search(r"SRID=(\d+);", wkt)
    if srid_match:
        srid = int(srid_match.group(1))
    else:
        raise ValueError("Input does not have a valid SRID declaration.")
    
    # Extract the content inside the parentheses
    content = re.search(r"LINESTRING\s*\((.*)\)", wkt).group(1)
    
    # Split by comma and then split each segment by space to get X, Y, Z
    points = [tuple(map(float, segment.split())) for segment in content.split(",")]
    
    start_point = Point(points[0], srid)
    end_point = Point(points[1], srid)
    
    linez = LineStringZ(start_point, end_point)
    
    return linez


def get_query(x1, y1, x2, y2, srid):
    with open("test_query.sql", "r") as f:
        query = f.read()
        query = query.replace("minX", str(x1))
        query = query.replace("minY", str(y1))
        query = query.replace("maxX", str(x2))
        query = query.replace("maxY", str(y2))
        query = query.replace("srid", str(srid))
        # Temporary
        query = query.replace("fkb_bygning", "fkb_bygning_42")

    return query


    """
        feature = dict(type="Feature", properties=dict(id=building_id), geometry=object)
        features.append(feature)
    
    return features
    """ 

def extract_objects(objects, building_id):
    features = list()
    
    for object in objects:
        object = json.loads(object)
        object.pop("crs", None)
        
        feature = dict(type="Feature", properties=dict(id=building_id), geometry=object)
        features.append(feature)
    
    return features

    
def extract_points(objects, building_id):
    points = list()
    
    
    for object in objects:
        object = json.loads(object)
        object.pop("crs", None)        
        
        # Go through all coordinates
        for coordinate in object["coordinates"]:
            # Round for duplicate removal later
            # coordinate = [round(c) for c in coordinate]
            points.append(coordinate)
        
    # Remove duplicates
    points = list(set([list(p) for p in points]))
    
    features = list()
    
    # Create point geometry, then add as feature with building id
    for point in points:
        point = dict(type="Point", coordinates=point)
        feature = dict(type="Feature", properties=dict(id=building_id), geometry=point)
        features.append(feature)
    return features


from collections import defaultdict


class Node:
    def __init__(self, coords):
        self.coords = coords
        self.edge_nodes = []
        
    def add_edge(self, node):
        if node not in self.edge_nodes:
            self.edge_nodes.append(node)
        if self not in node.edge_nodes:
            node.edge_nodes.append(self)
        
    def __eq__(self, other):
        return self.coords == other.coords

    def __hash__(self):
        return hash(tuple(self.coords))
        
    def __repr__(self):
        return f"({self.coords})"


class Polygon:
    def __init__(self, nodes):
        self.nodes = nodes
    
    def get_nodes(self):
        return [node.coords for node in self.nodes]
        
    def __eq__(self, other):
        # Check if two polygons have the same set of nodes.
        return set(self.nodes) == set(other.nodes)

    def __iter__(self):
        return iter(self.nodes)
        
    def __repr__(self):
        return f"Polygon({', '.join([str(node) for node in self.nodes])})"


class Graph:
    def __init__(self):
        self.nodes = []
        self.nodes_with_edges = {}
        self.used_edges = set()
        self.polygons = []

    def _get_or_create_node(self, coords):
        node = next((n for n in self.nodes if n.coords == coords), None)
        if node is None:
            node = Node(coords)
            self.nodes.append(node)
        return node

    def add_edge(self, coords1, coords2):
        node1 = self._get_or_create_node(coords1)
        node2 = self._get_or_create_node(coords2)
        node1.add_edge(node2)

    def retrieve_polygons(self):
        for start_node in self.nodes:
            visited = set()  # Reset visited nodes for each starting node
            stack = [(start_node, [start_node])]
            
            while stack:
                current, path = stack.pop()
                visited.add(current)
                
                for neighbor in current.edge_nodes:
                    if neighbor == path[0] and len(path) > 2:  # cycle found
                        new_polygon = Polygon(path + [path[0]])  # Close the loop by appending starting node
                        if new_polygon not in self.polygons:
                            self.polygons.append(new_polygon)
                    elif neighbor not in path and neighbor not in visited:
                        new_path = list(path)
                        new_path.append(neighbor)
                        stack.append((neighbor, new_path))

    def retrieve_polygons(self):
        for start_node in self.nodes:
            stack = [(start_node, [start_node])]
            
            while stack:
                current, path = stack.pop()
                
                for neighbor in current.edge_nodes:
                    if neighbor == path[0] and len(path) > 2:  # cycle found
                        new_polygon = Polygon(path + [path[0]])  # Close the loop by appending starting node
                        if new_polygon not in self.polygons:
                            self.polygons.append(new_polygon)
                    elif neighbor not in path:
                        new_path = list(path)
                        new_path.append(neighbor)
                        stack.append((neighbor, new_path))
    
    def retrieve_polygons(self):
        for start_node in self.nodes:
            for neighbor in start_node.edge_nodes:
                visited = set([start_node])
                stack = [(neighbor, [start_node, neighbor])]

                while stack:
                    current, path = stack.pop()
                    visited.add(current)

                    for next_node in current.edge_nodes:
                        if next_node == path[0] and len(path) > 2:  # cycle found
                            new_polygon = Polygon(path + [path[0]])  # Close the loop by appending starting node
                            if new_polygon not in self.polygons:
                                self.polygons.append(new_polygon)
                        elif next_node not in path and next_node not in visited:
                            new_path = list(path)
                            new_path.append(next_node)
                            stack.append((next_node, new_path))
    
    def __repr__(self):
        return f"Graph with {len(self.nodes)} nodes and {len(self.polygons)} polygons"


def canonical_form(polygon):
    """
    Convert the polygon into its canonical form.
    """
    coords_list = [tuple(node.coords) for node in polygon.nodes]
    min_index = min(range(len(coords_list)), key=lambda i: coords_list[i])
    rotated_coords = tuple(coords_list[min_index:] + coords_list[:min_index])
    
    return rotated_coords

def remove_duplicate_polygons(polygons):
    seen = set()
    unique_polygons = []

    for polygon in polygons:
        canon = canonical_form(polygon)
        if canon not in seen:
            seen.add(canon)
            unique_polygons.append(polygon)

    return unique_polygons

def is_subset(polygon_a, polygon_b):
    """
    Check if polygon_a is a subset of polygon_b.
    """
    return all(point in polygon_b for point in polygon_a)

def filter_polygons(polygons):
    """
    Filter out polygons that are supersets of another polygon.
    """
    filtered_polygons = []
    for polygon in polygons:
        if not any(is_subset(other, polygon) for other in polygons if other != polygon):
            filtered_polygons.append(polygon)
    return filtered_polygons

def polygons_to_geojson(polygons):
    """
    Convert a list of polygons to a GeoJSON feature collection.
    """
    geojson = {
        "type": "FeatureCollection",
        "features": []
    }

    for polygon in polygons:
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [polygon.get_nodes()]
            },
            "properties": {}
        }
        geojson['features'].append(feature)

    return geojson


if __name__ == "__main__":
    
    # Get the database cursor for retrieving data
    cur = get_cursor()

    # Create the query with the correct bounding box and srid
    query = get_query(482000, 6471000, 483000, 6472000, 25832)
    
    # Execute the query
    cur.execute(query)
    
    # Get all of the results
    results = cur.fetchall()
    
    # geojson for storing all of the data
    geojson = dict(type="FeatureCollection", features=[])
    
    polygons = []
    
    # The results are currently grouped per building, therefore we have a unique building
    # id for each building and we know which polygons correspond to which building
    for residx, result in enumerate(results):
        
        feature = dict(type="Feature", properties=dict(id=result[0]), geometry=[])
        point = dict(type="Point", coordinates=None)
        
        building_id = result[0]
        geometries = result[1:]
        
        building_omrade = geometries[0]
        monelinjer = geometries[1]
        takkanter = geometries[2]
        bygningslinjer = geometries[3]
        taksprang = geometries[4]
        
        
        # Parse the geometries
        #building_omrade = json.loads(building_omrade)
        # Remove crs from dictionary
        #building_omrade.pop("crs", None)
        # Add to geometries
        #geojson["features"].append(building_omrade)
        
        if monelinjer != None:
            #monelinje_points = extract_points(monelinjer, building_id)
            monelinje_objects = extract_objects(monelinjer, building_id)
            #geojson["features"].extend(monelinje_points)
            geojson["features"].extend(monelinje_objects)
        else:
            print("monelinjer is None")
        
        if takkanter != None:
            #takkant_points = extract_points(takkanter, building_id)
            takkant_objects = extract_objects(takkanter, building_id)
            #geojson["features"].extend(takkant_points)
            geojson["features"].extend(takkant_objects)
        else:
            print("takkanter is None")
        
        if bygningslinjer != None:
            #bygningslinje_points = extract_points(bygningslinjer, building_id)
            bygningslinje_objects = extract_objects(bygningslinjer, building_id)
            #geojson["features"].extend(bygningslinje_points)
            geojson["features"].extend(bygningslinje_objects)
        else:
            print("bygningslinjer is None")
        
        if taksprang != None:
            #takspring_points = extract_points(taksprang, building_id)
            takspring_objects = extract_objects(taksprang, building_id)
            #geojson["features"].extend(takspring_points)
            geojson["features"].extend(takspring_objects)
        else:
            print("taksprang is None")
            
        # data = geojson["features"]
        
        # geojson["features"] = []
        
        
        def round_coords(coords):
            return [round(c, 0) for c in coords]

        data = geojson["features"]
        geojson["features"] = []

        graph = Graph()
        
        for feature in data:
            coords = feature['geometry']['coordinates']
            for coordx in range(len(coords) - 1):
                graph.add_edge(coords[coordx], coords[coordx + 1])

        for node in graph.nodes:
            print("Node:", node, " | Edges:", [edge.coords for edge in node.edge_nodes])
        
        graph.retrieve_polygons()
        
        for p in graph.polygons:
            print(p, "\n")
        
        unique_polygons = filter_polygons(graph.polygons)
        
        print("\nUNIQUE\n")
        
        for p in unique_polygons:
            print(p, "\n")
        
        polygons.extend(unique_polygons)
        
    geojson_string = polygons_to_geojson(polygons)
    
    with open(f"linjer_{residx}.json", "w") as f:
        json.dump(geojson_string, f)