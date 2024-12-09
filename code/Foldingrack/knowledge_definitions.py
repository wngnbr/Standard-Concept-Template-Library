from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def hook_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Regular_hook)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = np.array([
                -obj.base_size[0] - obj.middle_size[0] * np.cos(obj.middle_rotation[0]) + obj.middle_size[1] / 2 * np.sin(obj.middle_rotation[0]),
                obj.middle_size[0] * np.sin(obj.middle_rotation[0]) + obj.middle_offset[0] - obj.circle_radius[0] + obj.middle_size[1] / 2 * np.cos(obj.middle_rotation[0]),
                obj.middle_offset[1]
            ])
            mesh_rotation = np.array([np.pi / 2, 0, -np.pi / 2])

            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            return (_pt_radius >= obj.circle_radius[0] - obj.middle_size[1] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.circle_radius[0] + AFFORDACE_PROXIMITY_THRES and
                    _pt[1] >= -obj.middle_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.middle_size[2]/2 + AFFORDACE_PROXIMITY_THRES)
        
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
