from concept_template import *
from geometry_template import *
from knowledge_utils import *


def cover_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Simplified_Cover)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = np.array([0, obj.size[1]/2, obj.size[2]/2])
            mesh_rotation = np.array([0, 0, 0])
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            return (_pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[2] >= -obj.size[2]/3 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Carved_Cover)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = np.array([0, (obj.outer_size[1] + obj.inner_size[1]) / 2, obj.outer_size[2] / 2])
            mesh_rotation = np.array([0, 0, 0])
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            return (_pt[0] >= -obj.outer_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[0] <= obj.outer_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -(obj.outer_size[1] - obj.inner_size[1])/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= (obj.outer_size[1] - obj.inner_size[1])/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[2] >= -obj.outer_size[2]/3 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[2] <= obj.outer_size[2]/2 + AFFORDACE_PROXIMITY_THRES)
    
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
