from concept_template import *
from geometry_template import *
from knowledge_utils import *


def lid_affordance(obj, pt):

    def is_affordance(obj, pt):
        if (isinstance(obj, Cylindrical_Lid)):
            _pt = inverse_transformation(pt, obj.position, obj.rotation)

            lid_y_offset = -(obj.outer_size[2] - obj.inner_size[2]) / 2
            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            return (_pt_radius >= obj.inner_size[1] - AFFORDACE_PROXIMITY_THRES and
                    _pt_radius <= obj.outer_size[0] + AFFORDACE_PROXIMITY_THRES and
                    _pt[1] >= -obj.inner_size[2]/2 + lid_y_offset - AFFORDACE_PROXIMITY_THRES and
                    _pt[1] <= obj.outer_size[2]/2 + AFFORDACE_PROXIMITY_THRES)
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
