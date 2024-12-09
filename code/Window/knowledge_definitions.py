from concept_template import *
from geometry_template import *
from knowledge_utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):
        if (isinstance(obj, Cuboidal_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]
            
            for i in range(obj.num_of_handle[0]):
                position_z = 0
                z_layer_position = obj.handle_z_position[i]
                if obj.window_type == 0:
                    if z_layer_position == -1:
                        position_z -= obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 0:
                        position_z += obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 1:
                        position_z += obj.windows_size["size_1"][1]
                elif obj.window_type == 1:
                    if z_layer_position >= 0:
                        position_z += obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 1:
                        position_z += obj.windows_size["size_1"][1]
                    if z_layer_position >= 2:
                        position_z += obj.windows_size["size_2"][1]
                    if z_layer_position == 3:
                        position_z += obj.windows_size["size_3"][1]
                handle_mesh_position = [obj.offset_x[i], 0, position_z + obj.size[2] / 2]
                handle_mesh_rotation = [0, 0, 0]
                __pt = inverse_transformation(_pt, handle_mesh_position, handle_mesh_rotation)
            
                if (__pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
                
            for i in range(obj.num_of_handle[1]):
                position_z = 0
                z_layer_position = obj.handle_z_position[i]
                if obj.window_type == 0:
                    if z_layer_position == -1:
                        position_z -= (
                            obj.windows_size["size_0"][1] / 2 + obj.windows_size["size_2"][1]
                        )
                    if z_layer_position >= 0:
                        position_z += obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 1:
                        position_z += obj.windows_size["size_0"][1]
                elif obj.window_type == 1:
                    if z_layer_position >= 0:
                        position_z -= obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 1:
                        position_z += obj.windows_size["size_0"][1]
                    if z_layer_position >= 2:
                        position_z += obj.windows_size["size_1"][1]
                    if z_layer_position == 3:
                        position_z += obj.windows_size["size_2"][1]
                handle_mesh_position = [obj.offset_x[i+2], 0, position_z - obj.size[2] / 2]
                handle_mesh_rotation = [0, 0, 0]
                __pt = inverse_transformation(_pt, handle_mesh_position, handle_mesh_rotation)
            
                if (__pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
                
            return False
        

        elif (isinstance(obj, Arched_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]
            
            central_angle = (np.arcsin((obj.seperation[0] / 2 + obj.bottom_size[1]) / obj.outer_size[0]) * 2)
            arch_offset_z = obj.bottom_size[2] - np.cos(central_angle / 2) * obj.outer_size[0]
            
            for i in range(obj.num_of_handle[0]):
                position_z = 0
                z_layer_position = obj.handle_z_position[i]
                if obj.window_type == 0:
                    if z_layer_position == -1:
                        position_z -= obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 0:
                        position_z += obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 1:
                        position_z += obj.windows_size["size_1"][1]
                elif obj.window_type == 1:
                    if z_layer_position >= 0:
                        position_z += obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 1:
                        position_z += obj.windows_size["size_1"][1]
                    if z_layer_position >= 2:
                        position_z += obj.windows_size["size_2"][1]
                    if z_layer_position == 3:
                        position_z += obj.windows_size["size_3"][1]
                main_mesh_position = np.array([
                    obj.offset_x[i],
                    0,
                    position_z + arch_offset_z,
                ])
                main_mesh_rotation = np.array([0, -np.pi / 2 + central_angle / 2, np.pi / 2])
                __pt = apply_transformation(np.array([_pt]), -main_mesh_position, -main_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
                _pt_radius = np.linalg.norm(__pt[[0, 2]])
                _pt_angle = np.arctan2(__pt[2], __pt[0])
                
                if (_pt_radius >= obj.seperation[0] / 2 / np.sin(central_angle / 2) + obj.thinner_handle[0] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.outer_size[0] + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.bottom_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.bottom_size[0]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
            
            for i in range(obj.num_of_handle[1]):
                position_z = 0
                z_layer_position = obj.handle_z_position[i]
                if obj.window_type == 0:
                    if z_layer_position == -1:
                        position_z -= (
                            obj.windows_size["size_0"][1] / 2 + obj.windows_size["size_2"][1]
                        )
                    if z_layer_position >= 0:
                        position_z += obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 1:
                        position_z += obj.windows_size["size_0"][1]
                elif obj.window_type == 1:
                    if z_layer_position >= 0:
                        position_z -= obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 1:
                        position_z += obj.windows_size["size_0"][1]
                    if z_layer_position >= 2:
                        position_z += obj.windows_size["size_1"][1]
                    if z_layer_position == 3:
                        position_z += obj.windows_size["size_2"][1]
                main_mesh_position = np.array([
                    obj.offset_x[i+2],
                    0,
                    position_z - arch_offset_z,
                ])
                main_mesh_rotation = np.array([0, np.pi / 2 + central_angle / 2, np.pi / 2])
                __pt = apply_transformation(np.array([_pt]), -main_mesh_position, -main_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
                _pt_radius = np.linalg.norm(__pt[[0, 2]])
                _pt_angle = np.arctan2(__pt[2], __pt[0])
                
                if (_pt_radius >= obj.seperation[0] / 2 / np.sin(central_angle / 2) + obj.thinner_handle[0] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.outer_size[0] + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.bottom_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.bottom_size[0]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
                
            return False
        

        elif (isinstance(obj, LShaped_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]
            
            
            for i in range(obj.num_of_handle[0]):
                position_z = 0
                z_layer_position = obj.handle_z_position[i]
                if obj.window_type == 0:
                    if z_layer_position == -1:
                        position_z -= obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 0:
                        position_z += obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 1:
                        position_z += obj.windows_size["size_1"][1]
                elif obj.window_type == 1:
                    if z_layer_position >= 0:
                        position_z += obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 1:
                        position_z += obj.windows_size["size_1"][1]
                    if z_layer_position >= 2:
                        position_z += obj.windows_size["size_2"][1]
                    if z_layer_position == 3:
                        position_z += obj.windows_size["size_3"][1]
                top_mesh_position = np.array([
                    obj.offset_x[i],
                    obj.offset_middle_y[0] + obj.offset_top_y[0],
                    position_z + obj.size_bottom[2] + obj.size_middle[2] + obj.size_top[2] / 2,
                ])
                top_mesh_rotation = np.array([0, 0, 0])
                __pt = apply_transformation(np.array([_pt]), -top_mesh_position, -top_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
                if (__pt[0] >= -obj.size_top[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[0] <= obj.size_top[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.size_top[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] >= -obj.size_top[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] <= obj.size_top[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
            
            for i in range(obj.num_of_handle[1]):
                position_z = 0
                z_layer_position = obj.handle_z_position[i]
                if obj.window_type == 0:
                    if z_layer_position == -1:
                        position_z -= (
                            obj.windows_size["size_0"][1] / 2 + obj.windows_size["size_2"][1]
                        )
                    if z_layer_position >= 0:
                        position_z += obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 1:
                        position_z += obj.windows_size["size_0"][1]
                elif obj.window_type == 1:
                    if z_layer_position >= 0:
                        position_z -= obj.windows_size["size_0"][1] / 2
                    if z_layer_position >= 1:
                        position_z += obj.windows_size["size_0"][1]
                    if z_layer_position >= 2:
                        position_z += obj.windows_size["size_1"][1]
                    if z_layer_position == 3:
                        position_z += obj.windows_size["size_2"][1]
                top_mesh_position = np.array([
                    obj.offset_x[i+2],
                    obj.offset_middle_y[0] + obj.offset_top_y[0],
                    position_z - obj.size_bottom[2] - obj.size_middle[2] - obj.size_top[2] / 2,
                ])
                top_mesh_rotation = np.array([0, 0, 0])
                __pt = apply_transformation(np.array([_pt]), -top_mesh_position, -top_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
                if (__pt[0] >= -obj.size_top[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[0] <= obj.size_top[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.size_top[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] >= -obj.size_top[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] <= obj.size_top[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
                
            return False
        

        return False
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
