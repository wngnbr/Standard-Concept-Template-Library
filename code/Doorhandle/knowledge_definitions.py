from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Regular_handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ')[0]

            main_mesh_position = np.array([
                obj.interpiece_offset_1[0] + obj.interpiece_offset_2[0],
                obj.interpiece_offset_1[1] + obj.interpiece_offset_2[1],
                obj.fixed_part_size[2] + obj.vertical_movable_size[2] + obj.horizontal_movable_size[2] / 2,
            ])
            main_mesh_rotation = np.array([0, 0, 0])
            __pt = inverse_transformation(_pt, main_mesh_position, main_mesh_rotation)

            return (__pt[0] >= -obj.horizontal_movable_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[0] <= obj.horizontal_movable_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= - AFFORDACE_PROXIMITY_THRES and
                __pt[1] <= obj.horizontal_movable_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                __pt[2] >= -obj.horizontal_movable_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[2] <= obj.horizontal_movable_size[2]/2 + AFFORDACE_PROXIMITY_THRES)
        


        elif (isinstance(obj, Knob_handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ')[0]

            main_mesh_position = np.array([0, 0, obj.fixed_part_size[1] + obj.sub_size[1] + obj.main_size[1] / 2])
            main_mesh_rotation = np.array([np.pi/2, 0, 0])
            __pt = inverse_transformation(_pt, main_mesh_position, main_mesh_rotation)

            _pt_radius = np.linalg.norm(__pt[[0, 2]])

            return (_pt_radius <= obj.main_size[0] + AFFORDACE_PROXIMITY_THRES and  
                    __pt[1] >= -obj.main_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.main_size[1]/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, TShaped_handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ')[0]

            main_mesh_position = np.array([
                obj.interpiece_offset[0],
                obj.interpiece_offset[1],
                obj.sub_size[1] + obj.main_size[2] / 2,
            ])
            main_mesh_rotation = np.array([0, 0, 0])
            __pt = inverse_transformation(_pt, main_mesh_position, main_mesh_rotation)

            return (__pt[0] >= -obj.main_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= obj.main_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.main_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.main_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    __pt[2] >= -obj.main_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= obj.main_size[2]/2 + AFFORDACE_PROXIMITY_THRES)
        
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
