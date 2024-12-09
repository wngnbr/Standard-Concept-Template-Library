from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def nozzle_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Cylindrical_cap)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            nozzle_mesh_position = np.array([0, obj.inner_size[2] / 2, 0])
            nozzle_mesh_rotation = np.array([0, 0, 0])
            __pt = inverse_transformation(_pt, nozzle_mesh_position, nozzle_mesh_rotation)

            _pt_radius = np.linalg.norm(__pt[[0, 2]])

            return (_pt_radius <= obj.outer_size[0] and
                    __pt[1] >= - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= (obj.outer_size[2] - obj.inner_size[2])/2 + AFFORDACE_PROXIMITY_THRES)
        
        

        elif (isinstance(obj, Regular_cap)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            _pt_radius = np.linalg.norm(__pt[[0, 2]])

            return (_pt_radius <= obj.size[0] and
                    __pt[1] >= - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES)
    
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
