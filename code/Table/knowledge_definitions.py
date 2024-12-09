from concept_template import *
from geometry_template import *
from knowledge_utils import *


def door_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Regular_door)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            for door_idx in range(obj.number_of_door[0]):
                mesh_position = np.array([
                    obj.door_offset[door_idx][0] + obj.handle_offset[door_idx][0] * np.cos(obj.door_rotation[door_idx]) + obj.handle_size[door_idx][2] / 2 * np.sin(obj.door_rotation[door_idx]),
                    obj.door_offset[door_idx][1] - obj.door_size[door_idx][1] / 2 + obj.handle_offset[door_idx][1],
                    obj.door_offset[door_idx][2] - obj.handle_offset[door_idx][0] * np.sin(obj.door_rotation[door_idx]) + obj.handle_size[door_idx][2] / 2 * np.cos(obj.door_rotation[door_idx])
                ])
                mesh_rotation = np.array([0, obj.door_rotation[door_idx], obj.handle_rotation[door_idx]])
                __pt = apply_transformation(np.array([_pt]), -mesh_position, -mesh_rotation, rotation_order='ZYX', offset_first=True)[0]

                if (__pt[0] >= -obj.handle_size[door_idx][0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[0] <= obj.handle_size[door_idx][0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.handle_size[door_idx][1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.handle_size[door_idx][1]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] >= -obj.handle_size[door_idx][2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] <= obj.handle_size[door_idx][2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
            
            return False
    
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def drawer_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Regular_drawer)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            for drawer_idx in range(obj.number_of_drawer[0]):
                for mesh_idx in range(6, 6 + obj.number_of_handle[drawer_idx]):
                    if obj.number_of_handle[drawer_idx] == 2:
                        position_sign = 1 if mesh_idx == 6 else -1
                    else:
                        position_sign = 0
                    mesh_position = [
                        obj.drawer_offset[drawer_idx][0] + obj.handle_offset[drawer_idx][0] + position_sign * obj.handle_separation[drawer_idx] / 2,
                        obj.drawer_offset[drawer_idx][1] + obj.handle_offset[drawer_idx][1],
                        obj.drawer_offset[drawer_idx][2] + obj.drawer_size[drawer_idx][2] / 2 + obj.front_size[drawer_idx][2] + obj.front_size[drawer_idx][2] / 2
                    ]
                    mesh_rotation = [0, 0, 0]
                    __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

                    if (__pt[0] >= -obj.handle_sizes[drawer_idx][0]/2 - AFFORDACE_PROXIMITY_THRES and 
                        __pt[0] <= obj.handle_sizes[drawer_idx][0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= -obj.handle_sizes[drawer_idx][1]/2 - AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] <= obj.handle_sizes[drawer_idx][1]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[2] >= -obj.handle_sizes[drawer_idx][2]/2 - AFFORDACE_PROXIMITY_THRES and 
                        __pt[2] <= obj.handle_sizes[drawer_idx][2]/2 + AFFORDACE_PROXIMITY_THRES):
                        return True
            
            return False
    
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
