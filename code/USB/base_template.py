import numpy as np
from knowledge_utils import *

class GeometryTemplate:

    def __init__(self, position, rotation, rotation_order):
        self.position = position
        self.rotation = rotation
        self.rotation_order = rotation_order


class ConceptTemplate:

    def __init__(self, position, rotation):
        self.position = position
        self.rotation = rotation

    def proximation(self, pt):
        dist = np.linalg.norm(self.overall_obj_pts - pt, axis=1)
        return np.min(dist) < PROXIMITY_THRES