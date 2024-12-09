from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def screen_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Regular_Screen)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY', offset_first=True)[0]

            mesh_position = np.array([
                0, 
                obj.offset[0] + obj.size[1] * np.cos(obj.screen_rotation[0]) / 2,
                obj.offset[1]
            ])
            mesh_rotation = np.array([0, 0, 0])
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            return (_pt[0] >= -obj.size[0]/4 - AFFORDACE_PROXIMITY_THRES and
                    _pt[0] <= obj.size[0]/4 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= obj.size[1]*3/8 - AFFORDACE_PROXIMITY_THRES and
                    _pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    _pt[2] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[2] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES)

        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
