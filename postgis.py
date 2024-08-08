from postgisutils import get_cursor, close_cursor
from graph import Graph, Node
import json

from tqdm import tqdm
import argparse

from utils import generate_query_boxes
from new_utils import subtract_smaller_polygons_from_larger, extract_unique_edges, build_graph, find_unique_polygons, convert_to_geojson_feature, validate_polygons, filter_redundant_polygons, merge_close_nodes, filter_duplicate_polygons

def is_collinear(node1, node2, node3, epsilon=1e-6):
    x1, y1 = node1.x, node1.y
    x2, y2 = node2.x, node2.y
    x3, y3 = node3.x, node3.y
    determinant = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
    return abs(determinant) < epsilon

def simplify_graph(nodes):
    """
    Simplify a graph by removing nodes that are collinear with their two neighbors,
    except those nodes that are critical junctions (having more or less than two neighbors).
    """
    count = 0
    removable_nodes = [node for node in nodes if node.is_collinear_with_neighbors()]
    for node in removable_nodes:
        count += 1
        # Remove the node by connecting its neighbors to each other
        n1, n2 = list(node.neighbors)
        n1.neighbors.remove(node)
        n2.neighbors.remove(node)
        n1.neighbors.add(n2)
        n2.neighbors.add(n1)
        nodes.remove(node)

    print("Simplified {} nodes".format(count))

    return nodes

def get_query(x1, y1, x2, y2, srid, box_srid):
    with open("query.sql", "r") as f:
        query = f.read()
        query = query.replace("minX", str(x1))
        query = query.replace("minY", str(y1))
        query = query.replace("maxX", str(x2))
        query = query.replace("maxY", str(y2))
        query = query.replace("box_srid", str(box_srid))
        query = query.replace("srid", str(srid))
        # Temporary
        # query = query.replace("fkb_bygning", "fkb_bygning_42") # Only choose Agder
        query = query.replace("tolerance_distance", str(0.01))

    return query

def data_to_geojson(polygons):
    """Generate a GeoJSON object for a list of polygons."""
    features = polygons
    return {
        "type": "FeatureCollection",
        "features": features
    }

def clean_geometry(geom):
    if not geom.is_valid:
        clean_geom = geom.buffer(0)
        if clean_geom.is_valid:
            return clean_geom
    return geom


def remove_area_based_redundant_polygons(polygons):
    # Sort polygons by area in descending order to start with the largest polygons
    polygons_sorted = sorted(polygons, key=lambda x: x['area'], reverse=True)

    removed_indices = set()
    for i, larger_polygon in enumerate(polygons_sorted):
        if i in removed_indices:
            continue  # Skip already marked polygons
        
        possible_combinations = [p for j, p in enumerate(polygons_sorted) if j != i]
        for j, smaller_polygon in enumerate(possible_combinations):
            if larger_polygon['area'] == sum(p['area'] for p in possible_combinations if p != smaller_polygon):
                # Mark the larger polygon for removal if its area is exactly matched by the combination of smaller polygons
                removed_indices.add(i)
                break

    # Create a new list excluding the redundant (larger) polygons
    remaining_polygons = [p for i, p in enumerate(polygons_sorted) if i not in removed_indices]

    return remaining_polygons
def extract_data_from_query(results):

    geojson_features = list()
    node_neighbor_dicts = list()
    training_query_polygons = list()

    already_seen_buildings = list()

    for residx, result in enumerate(results):
        
        building_id = result[0]
        building_box = result[1]

        if building_id in already_seen_buildings:
            continue
        
        already_seen_buildings.append(building_id)

        json_geometries = list()

        for value in result[2].values():
            json_geometries.append(json.loads(value))

        building_box = json.loads(building_box)

        json_geometries = list(set([json.dumps(jg) for jg in json_geometries]))
        json_geometries = [json.loads(jg) for jg in json_geometries]

        crs = "EPSG:25833"

        margin = 0.001

        # Example usage:
        # edges = extract_edges(json_geometries)

        edges = extract_unique_edges(json_geometries)

        # edges = merge_collinear_edges(edges)

        # edges = merge_close_nodes(edges, tolerance=margin)

        graph = build_graph(edges)

        polygons = find_unique_polygons(graph)

        polygons = validate_polygons(polygons)

        # polygons = filter_duplicate_polygons(polygons)

        # polygons = filter_redundant_polygons(polygons)

        polygons = subtract_smaller_polygons_from_larger(polygons)


        geojson_feature = [convert_to_geojson_feature(polygon) for polygon in polygons]

        geojson_features.extend(geojson_feature)
        

    return node_neighbor_dicts, geojson_features


def parse_args():
    parser = argparse.ArgumentParser(description='PostGIS Extract')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    return parser.parse_args()

if __name__ == "__main__":    
    
    # Get the database cursor for retrieving data
    cur = get_cursor()
    
    args = parse_args()
    
    bounding_box = (650352, 7728709, 652552, 7729909)
    
    # bounding_box = (443966.4976,6444679.4794,445589.1116,6447090.2756)

    # Ø 650801 - Ø 655857
    #N 7728992 - N 7739380

    # bounding_box = (650801, 7728992, 655857, 7739380)

    # Split the bounding box into smaller chunks to enable faster queries
    # The smaller the chunks, the faster the queries, but the more queries we have to make
    # The larger the chunks, the slower the queries, but the less queries we have to make
    # The optimal chunk size is dependent on the data and the database
    # For this dataset, a chunk size of 100x100 meters seems to be optimal
    box_size = [100, 100]
    
    boxes = generate_query_boxes(*bounding_box, *box_size)
    
    if args.verbose:
        print("Bounding box: {}".format(bounding_box))
        print("Box size: {}".format(box_size))
        print("Number of boxes: {}".format(len(boxes)))
    
    # Create a list of all of the queries we need to make
    geojson_polygons = list()
    training_geojson_polygons = list()
    
    for boxid, box in tqdm(enumerate(boxes), total=len(boxes)):
        query = get_query(*box, srid=4326, box_srid=25833)
        
        if args.verbose:
            print("Running query for chunk ({}, {})".format(box[0], box[3]))
        
        try:
            cur.execute(query)
        except Exception as e:
            cur = close_cursor(cur)
            cur = get_cursor()
            print("Error in query: {}, for bbox {}, srid {}, and box_srid {}".format(e, box, 4326, 25833))
            continue
        
        result = cur.fetchall()
        
        # print("Query took {} seconds".format(end - start))
        
        if args.verbose:
            print("Found {} buildings - Extracting data".format(len(result)))
        
        node_neighbors, geojson = extract_data_from_query(result)

        if len(geojson) == 0:
            continue

        if args.verbose:
            print("Finished extracting data")
        
        geojson_polygons.extend(geojson)

    geojson_polygons = data_to_geojson(geojson_polygons)

    print(geojson_polygons)

    with open("geojson.json", "w") as f1:
        f1.write(json.dumps(geojson_polygons))

    close_cursor(cur)


    # TODO: FIX