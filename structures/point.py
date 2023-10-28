
class Point:
    def __init__(self, data: tuple, srid: int, building_id) -> None:
        self.x = data[0]
        self.y = data[1]
        self.z = data[2]
        self.srid = srid
        self.building_id = building_id
        
        self.nodes = list()
    
    def __str__(self) -> str:
        return f"({round(self.x, 1)}, {round(self.y, 1)}, {round(self.z, 1)} ({self.srid}) ({self.building_id}))"