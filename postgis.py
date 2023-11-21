from postgisutils import get_cursor
from data_structs import Graph, Node, Edge, Polygon

import json

def get_query(x1, y1, x2, y2, srid):
    with open("test_query.sql", "r") as f:
        query = f.read()
        query = query.replace("minX", str(x1))
        query = query.replace("minY", str(y1))
        query = query.replace("maxX", str(x2))
        query = query.replace("maxY", str(y2))
        query = query.replace("srid", str(srid))
        # Temporary
        query = query.replace("fkb_bygning", "fkb_bygning_42") # Only choose Agder

    return query


def generate_geojson(data, type):
    """Generate a GeoJSON representation for a given polygon."""
    
    if type == "Polygon":
        coordinates = [[node.x, node.y, node.z] for node in data.nodes]
        coordinates.append(coordinates[0])  # Close the loop
        coordinates = [coordinates]
    elif type == "Point":
        coordinates = data.x, data.y, data.z
    else:
        raise ValueError("Invalid type: {}, choose from {}".format(type, ["Polygon", "Point"]))
    
    
    return {
        "type": "Feature",
        "geometry": {
            "type": f"{type}",
            "coordinates": coordinates
        },
        "properties": {}
    }

# Split edges that have a node directly on them
def split_edges(graph):
    edges_to_remove = []
    edges_to_add = []

    for edge in graph.edges:
        for node in graph.nodes:
            if node not in [edge.node1, edge.node2] and is_point_on_line(node, edge):
                edges_to_remove.append(edge)
                edges_to_add.append(Edge(edge.node1, node))
                edges_to_add.append(Edge(node, edge.node2))
                break  # Once we split an edge for one node, we don't want to check other nodes for the same edge

    # Update the graph's edges
    for edge in edges_to_remove:
        graph.edges.remove(edge)
    for edge in edges_to_add:
        graph.edges.append(edge)

    return graph
    
def is_point_on_line(node, edge):
    # This is a simple method to determine if a point lies on a line segment.
    # It checks if the sum of distances from the point to the line segment's endpoints is equal to the distance between the endpoints.
    d1 = distance(node, edge.node1)
    d2 = distance(node, edge.node2)
    d_total = distance(edge.node1, edge.node2)

    return abs(d1 + d2 - d_total) < 1e-2  # We use a small threshold to account for floating point inaccuracies.

def distance(node1, node2):
    return ((node1.x - node2.x)**2 + (node1.y - node2.y)**2 + (node1.z - node2.z)**2) ** 0.5


def polygons_to_geojson(polygons):
    """Generate a GeoJSON object for a list of polygons."""
    features = polygons
    return {
        "type": "FeatureCollection",
        "features": features
    }

def remove_duplicate_polygons(polygons):
    seen = set()
    unique_polygons = []

    for polygon in polygons:
        # Create a canonical representation by sorting nodes in the polygon
        # Using tuple of tuples as a representation to make it hashable
        canonical = tuple(sorted((node.x, node.y, node.z) for node in polygon))
        if canonical not in seen:
            seen.add(canonical)
            unique_polygons.append(polygon)
        else:
            print("Found duplicate polygon")

    return unique_polygons

def filter_polygons_by_area(polygons, min_area=1e-2):
    """Remove polygons with an area smaller than min_area."""
    
    polygons = [Polygon(polygon) for polygon in polygons]
    
    for p in polygons:
        p.calculate_area()
    
    num_before_filter = len(polygons)
    
    polygons = [polygon for polygon in polygons if polygon.area >= min_area]
    
    num_after_filter = len(polygons)
    
    if num_before_filter != num_after_filter:
        print("Filtered {} polygons with area smaller than {}".format(num_before_filter - num_after_filter, min_area))
    
    polygons = sorted(polygons, key=lambda p: p.area, reverse=True)

    return polygons

def share_edges(poly1, poly2):
    """Check if two polygons share any edges."""
    for i in range(len(poly1) - 1):
        edge1 = (poly1[i], poly1[i+1])
        for j in range(len(poly2) - 1):
            edge2 = (poly2[j], poly2[j+1])
            if edge1 == edge2 or edge1 == (edge2[1], edge2[0]):
                return True
    return False

def is_contained(smaller, larger):
    """Check if the smaller polygon is partly or wholly contained in the larger polygon."""
    # Note: For a more precise check, consider using a point-in-polygon algorithm.
    # For simplicity, here we're checking if any vertex of the smaller polygon is inside the larger polygon.
    for vertex in smaller:
        if point_in_polygon(vertex, larger):
            return True
    return False

def point_in_polygon(point, polygon):
    """Determine if a point is inside a polygon."""
    n = len(polygon)
    odd_nodes = False
    j = n - 1  # The last vertex is the previous one to compare with the first.

    for i in range(n):
        if polygon[i].y < point.y and polygon[j].y >= point.y or polygon[j].y < point.y and polygon[i].y >= point.y:
            if polygon[i].x + (point.y - polygon[i].y) / (polygon[j].y - polygon[i].y) * (polygon[j].x - polygon[i].x) < point.x:
                odd_nodes = not odd_nodes
        j = i

    return odd_nodes

def compute_area(polygon):
    """Compute the area of a 2D polygon using the Shoelace formula."""
    n = len(polygon) - 1  # The last vertex is the same as the first
    area = 0
    for i in range(n):
        j = (i + 1) % n
        area += polygon[i].x * polygon[j].y
        area -= polygon[j].x * polygon[i].y
    area = abs(area) / 2.0
    return area

def polygons_share_edge(p1, p2):
    """Check if two polygons share an edge."""
    for i in range(len(p1)-1):
        for j in range(len(p2)-1):
            if (p1[i] == p2[j] and p1[i+1] == p2[j+1]) or (p1[i] == p2[j+1] and p1[i+1] == p2[j]):
                return True
    return False

def is_fully_contained(inside, outside):
    """Check if every vertex of the first polygon is inside the second one."""
    return all([point_in_polygon(v, outside) for v in inside])

def filter_polygons(polygons):
    """Filter polygons based on given criteria."""
    to_remove = set()

    for i, poly1 in enumerate(polygons):
        for j, poly2 in enumerate(polygons):
            if i != j:
                if polygons_share_edge(poly1, poly2):
                    if compute_area(poly1) < compute_area(poly2) and is_fully_contained(poly1, poly2):
                        to_remove.add(tuple(poly2))
                    elif compute_area(poly2) < compute_area(poly1) and is_fully_contained(poly2, poly1):
                        to_remove.add(tuple(poly1))

    return [polygon for polygon in polygons if tuple(polygon) not in to_remove]

def extract_data_from_query(results):
    
    query_polygons = list()
    query_points = list()
    
    for residx, result in enumerate(results):
        
        building_id = result[0]
        geometries = result[1:]        
        
        graph = Graph()
        
        for gidx, geometry in enumerate(geometries):
            
            strings = list()
            
            if isinstance(geometry, str):
                strings.append(geometry)
            elif isinstance(geometry, list) and geometry[0] is not None:
                for g in geometry:
                    strings.append(g)
            else:
                continue
                
            for geometry_string in strings:
                
                geom = json.loads(geometry_string)

                nodes = geom["coordinates"]
                
                for nodex in range(len(nodes) - 1):
                    curr_node = Node(*nodes[nodex])
                    next_node = Node(*nodes[nodex + 1])
                    graph.add_node(curr_node)
                    graph.add_node(next_node)
                    graph.connect_nodes(curr_node, next_node)
        
        graph.remove_duplicate_nodes() # There exists a lot of duplicate nodes from the linestrings
        
        
        graph = split_edges(graph) # If a node is on an edge, then merge the node into the edge
        
        graph.populate_node_neighbours() # It is easier to traverse the graph using neighbor attribute, so populate it
        
        polygons = graph.find_polygons() # Find polygons in the graph
        
        polygons = remove_duplicate_polygons(polygons) # Remove all of the duplicate polygons (if any)
        
        polygons = filter_polygons_by_area(polygons, 1) # Remove all polygons smaller than a certain area
        
        #print("")
        
        #print("Found {} polygons for building {}".format(len(polygons), building_id))
        
        polygons = [generate_geojson(polygon, "Polygon") for polygon in polygons]
        points = [generate_geojson(node, "Point") for node in graph.nodes]
        
        query_polygons.extend(polygons)
        query_points.extend(points)
                
    return query_polygons, query_points

if __name__ == "__main__":
    
    # Get the database cursor for retrieving data
    cur = get_cursor()
    
    bounding_box = (481000, 6471000, 481100, 6471100)
    
    # Split the bounding box into smaller chunks to enable faster queries
    # The smaller the chunks, the faster the queries, but the more queries we have to make
    # The larger the chunks, the slower the queries, but the less queries we have to make
    # The optimal chunk size is dependent on the data and the database
    # For this dataset, a chunk size of 100x100 meters seems to be optimal
    chunk_size = 50
    
    # Create a list of all of the queries we need to make
    geojson_polygons = list()
    geojson_points = list()
    
    for x in range(bounding_box[0], bounding_box[2], chunk_size):
        for y in range(bounding_box[1], bounding_box[3], chunk_size):
            query = get_query(x, y, x + chunk_size, y + chunk_size, 25832)
            print("Running query for chunk ({}, {})".format(x, y))
            cur.execute(query)
            result = cur.fetchall()
            print("Found {} buildings - Extracting data".format(len(result)))
                
            polygons, points = extract_data_from_query(result)
                
            print("Finished extracting data")
                
            if polygons is None:
                continue
            if points is None:
                continue
            
            geojson_polygons.extend(polygons)
            geojson_points.extend(points)
            
    
    geojson_polygons = polygons_to_geojson(geojson_polygons)
    geojson_points = polygons_to_geojson(geojson_points)
    
    with open("polygons.json", "w") as f:
        f.write(json.dumps(geojson_polygons))
    
    with open("points.json", "w") as f:
        f.write(json.dumps(geojson_points))
    
    exit("")
    
    for query in queries:
        cur.execute(query)
        geojson_polygons
    
    # Create the query with the correct bounding box and srid
    query = get_query(481000, 6471000, 482000, 6472000, 25832)
    # Execute the query
    cur.execute(query)
    # Get all of the results
    results = cur.fetchall()
    
    # The results are currently grouped per building, therefore we have a unique building
    # id for each building and we know which polygons correspond to which building
    
    
    geojson_polygons = list()
    
    for residx, result in enumerate(results):
        
        building_id = result[0]
        geometries = result[1:]        
        
        graph = Graph()
        
        for gidx, geometry in enumerate(geometries):
            
            strings = list()
            
            if isinstance(geometry, str):
                strings.append(geometry)
            elif isinstance(geometry, list) and geometry[0] is not None:
                for g in geometry:
                    strings.append(g)
            else:
                continue
                
            for geometry_string in strings:
                
                geom = json.loads(geometry_string)

                nodes = geom["coordinates"]
                
                for nodex in range(len(nodes) - 1):
                    curr_node = Node(*nodes[nodex])
                    next_node = Node(*nodes[nodex + 1])
                    graph.add_node(curr_node)
                    graph.add_node(next_node)
                    graph.connect_nodes(curr_node, next_node)
        
        graph.remove_duplicate_nodes() # There exists a lot of duplicate nodes from the linestrings
        
        graph = split_edges(graph) # If a node is on an edge, then merge the node into the edge
        
        graph.populate_node_neighbours() # It is easier to traverse the graph using neighbor attribute, so populate it
        
        polygons = graph.find_polygons() # Find polygons in the graph
        
        polygons = remove_duplicate_polygons(polygons) # Remove all of the duplicate polygons (if any)
        
        polygons = filter_polygons_by_area(polygons, 1) # Remove all polygons smaller than a certain area
        
        print("")
        
        print("Found {} polygons for building {}".format(len(polygons), building_id))
        
        for polygon in polygons:
            geojson_polygons.append(generate_geojson(polygon))
            
    with open("polygons_2.json", "w") as f:
        f.write(json.dumps(polygons_to_geojson(geojson_polygons)))