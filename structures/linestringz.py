from structures.point import Point

from structures.utils import calc_3d_length

class LineStringZ:
    
    def __init__(self, start_point, end_point) -> None:
        self.start = start_point
        self.end = end_point
        self.length = self.get_length()
        
        assert self.start.srid == self.end.srid, "SRID of start and end point must be the same."
        self.srid = self.start.srid
    
    def get_length(self):
        return round(calc_3d_length(self.start, self.end), 4)