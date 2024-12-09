from concept_template import *
from geometry_template import *
from knowledge_utils import *


def lever_affordance(obj, pt):

    def is_affordance(obj, pt):
        if (isinstance(obj, Regular_lever)):
            # _pt = inverse_transformation(pt, obj.position, obj.rotation)
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), 'YXZ')[0]


            # left-handle
            handle_position = np.array([
                0, 
                obj.level_handle_offset[0] * np.cos(obj.level_handle_rotation[0]),
                (obj.level_support_size[2] + obj.level_handle_size[2]) / 2 - obj.level_handle_offset[0] * np.sin(obj.level_handle_rotation[0])
            ])
            handle_rotation = np.array([-obj.level_handle_rotation[0], 0, 0])
            __pt = inverse_transformation(_pt, handle_position, handle_rotation)

            if (__pt[0] >= -obj.level_handle_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[0] <= obj.level_handle_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.level_handle_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[1] <= -obj.level_handle_size[1]/4 + AFFORDACE_PROXIMITY_THRES and
                __pt[2] >= -obj.level_handle_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[2] <= obj.level_handle_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            

            # right-handle
            handle_position = np.array([
                0, 
                obj.level_handle_offset[0] * np.cos(obj.level_handle_rotation[0]),
                -(obj.level_support_size[2] + obj.level_handle_size[2]) / 2 + obj.level_handle_offset[0] * np.sin(obj.level_handle_rotation[0])
            ])
            handle_rotation = np.array([obj.level_handle_rotation[0], 0, 0])
            __pt = inverse_transformation(_pt, handle_position, handle_rotation)

            if (__pt[0] >= -obj.level_handle_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[0] <= obj.level_handle_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.level_handle_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[1] <= -obj.level_handle_size[1]/4 + AFFORDACE_PROXIMITY_THRES and
                __pt[2] >= -obj.level_handle_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[2] <= obj.level_handle_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
                
            return False

        return False
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
