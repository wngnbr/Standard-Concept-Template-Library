from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def cover_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Domed_Cover)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            bottom_mesh_position = np.array([
                0,
                obj.bottom_size[2] / 2,
                0,
            ])
            bottom_mesh_rotation = np.array([0, 0, 0])
            _pt = inverse_transformation(_pt, bottom_mesh_position, bottom_mesh_rotation)

            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            return (_pt[1] >= -obj.bottom_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.bottom_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= (obj.bottom_size[0] + obj.bottom_size[1]) / 2 - _pt[1] * (obj.bottom_size[1] - obj.bottom_size[0]) / obj.bottom_size[2] + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Cylindrical_Cover)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            radius_middle = obj.outer_size[0] + (obj.outer_size[2] - obj.inner_size[2]) * (obj.outer_size[1] - obj.outer_size[0]) / obj.outer_size[2]
            bottom_mesh_position = np.array([
                0,
                (obj.inner_size[2] + obj.outer_size[2]) / 2, 
                0,
            ])
            bottom_mesh_rotation = np.array([0, 0, 0])
            _pt = inverse_transformation(_pt, bottom_mesh_position, bottom_mesh_rotation)

            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            return (_pt[1] >= -(obj.outer_size[2] - obj.inner_size[2])/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= (obj.outer_size[2] - obj.inner_size[2])/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= (obj.outer_size[0] + radius_middle) / 2 - _pt[1] * (obj.outer_size[0] - radius_middle) / (obj.outer_size[2] - obj.inner_size[2]) + AFFORDACE_PROXIMITY_THRES)
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
