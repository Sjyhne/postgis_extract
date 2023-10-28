from structures.point import Point

import math

def calc_3d_length(start: Point, end: Point):
    
    return math.sqrt((end.x - start.x)**2 + (end.y - start.y)**2 + (end.z - start.z)**2)