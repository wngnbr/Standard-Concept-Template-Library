from concept_template import *
from geometry_template import *
from knowledge_utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Straight_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            behind_1_mesh_position_1 = np.array([-obj.behind_size[0] / 2, 0, -obj.behind_size[2] / 2])
            behind_1_mesh_rotation = np.array([0, -obj.handle_rotation[1], 0])
            behind_1_mesh_position_1 = np.array(adjust_position_from_rotation(behind_1_mesh_position_1, behind_1_mesh_rotation))
            behind_1_mesh_position_2 = np.array([
                obj.handle_separation[0] / 2 + obj.front_size[0] / 2 * np.cos(obj.handle_rotation[0]) + obj.front_size[2] * np.sin(obj.handle_rotation[0]) + obj.front_behind_offset[0], 
                obj.front_behind_offset[1],
                -obj.front_size[2] * np.cos(obj.handle_rotation[0]) + obj.front_size[0] / 2 * np.sin(obj.handle_rotation[0])
            ])
            behind_1_mesh_position = behind_1_mesh_position_1 + behind_1_mesh_position_2

            __pt = inverse_transformation(_pt, behind_1_mesh_position, behind_1_mesh_rotation)

            if (__pt[0] >= -obj.behind_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= obj.behind_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.behind_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.behind_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.behind_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.behind_size[2]/3 + AFFORDACE_PROXIMITY_THRES):
                return True
            

            behind_2_mesh_position_1 = np.array([obj.behind_size[0] / 2, 0, -obj.behind_size[2] / 2])
            behind_2_mesh_rotation = np.array([0, obj.handle_rotation[1], 0])
            behind_2_mesh_position_1 = np.array(adjust_position_from_rotation(behind_2_mesh_position_1, behind_2_mesh_rotation))
            behind_2_mesh_position_2 = np.array([
                -obj.handle_separation[0] / 2 - obj.front_size[0] / 2 * np.cos(obj.handle_rotation[0]) - obj.front_size[2] * np.sin(obj.handle_rotation[0]) - obj.front_behind_offset[0], 
                obj.front_behind_offset[1],
                -obj.front_size[2] * np.cos(obj.handle_rotation[0]) + obj.front_size[0] / 2 * np.sin(obj.handle_rotation[0])
            ])
            behind_2_mesh_position = behind_2_mesh_position_1 + behind_2_mesh_position_2

            __pt = inverse_transformation(_pt, behind_2_mesh_position, behind_2_mesh_rotation)

            if (__pt[0] >= -obj.behind_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= obj.behind_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.behind_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.behind_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.behind_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.behind_size[2]/3 + AFFORDACE_PROXIMITY_THRES):
                return True
            
            return False
    
                

        elif (isinstance(obj, Rear_Curved_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY', offset_first=True)[0]

            behind_1_mesh_rotation_1 = np.array([0, np.pi, np.pi])
            behind_1_mesh_position_1 = np.array([-obj.behind_size[1], 0, 0])
            behind_1_mesh_rotation_2 = np.array([0, -obj.handle_rotation[1], 0])
            behind_1_mesh_rotation_2_reverse = np.array([0, obj.handle_rotation[1], 0])
            behind_1_mesh_position_1 = np.array(adjust_position_from_rotation(behind_1_mesh_position_1, behind_1_mesh_rotation_2))
            behind_1_mesh_position_2 = np.array([
                obj.handle_separation[0] / 2 - obj.front_size[0] / 2 * np.cos(obj.handle_rotation[0]) + obj.front_size[2] * np.sin(obj.handle_rotation[0]) + obj.front_behind_offset[0], 
                obj.front_behind_offset[1],
                -obj.front_size[2] * np.cos(obj.handle_rotation[0]) + obj.front_size[0] / 2 * np.sin(obj.handle_rotation[0])
            ])
            behind_1_mesh_position = behind_1_mesh_position_1 + behind_1_mesh_position_2
            behind_1_mesh_rotation = behind_1_mesh_rotation_1 + behind_1_mesh_rotation_2_reverse

            __pt = apply_transformation(np.array([_pt]), -behind_1_mesh_position, -behind_1_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
            _pt_radius = np.linalg.norm(__pt[[0,2]])
            _pt_angle = np.arctan2(__pt[2], __pt[0])
            
            if (_pt_radius >= obj.behind_size[1] - AFFORDACE_PROXIMITY_THRES and 
                _pt_radius <= obj.behind_size[0] + AFFORDACE_PROXIMITY_THRES and
                __pt[1] >= -obj.behind_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.behind_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                _pt_angle >= obj.exist_angle[0]/8 and _pt_angle <= obj.exist_angle[0]):
                return True
        

            behind_2_mesh_rotation_1 = np.array([0, np.pi, 0])
            behind_2_mesh_position_1 = np.array([obj.behind_size[1], 0, 0])
            behind_2_mesh_rotation_2 = np.array([0, obj.handle_rotation[1], 0])
            behind_2_mesh_rotation_2_reverse = np.array([0, obj.handle_rotation[1], 0])
            behind_2_mesh_position_1 = np.array(adjust_position_from_rotation(behind_2_mesh_position_1, behind_2_mesh_rotation_2))
            behind_2_mesh_position_2 = np.array([
                -obj.handle_separation[0] / 2 + obj.front_size[0] / 2 * np.cos(obj.handle_rotation[0]) - obj.front_size[2] * np.sin(obj.handle_rotation[0]) - obj.front_behind_offset[0], 
                obj.front_behind_offset[1],
                -obj.front_size[2] * np.cos(obj.handle_rotation[0]) + obj.front_size[0] / 2 * np.sin(obj.handle_rotation[0])
            ])
            behind_2_mesh_position = behind_2_mesh_position_1 + behind_2_mesh_position_2
            behind_2_mesh_rotation = behind_2_mesh_rotation_1 + behind_2_mesh_rotation_2_reverse

            __pt = apply_transformation(np.array([_pt]), -behind_2_mesh_position, -behind_2_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
            _pt_radius = np.linalg.norm(__pt[[0,2]])
            _pt_angle = np.arctan2(__pt[2], __pt[0])
            
            if (_pt_radius >= obj.behind_size[1] - AFFORDACE_PROXIMITY_THRES and 
                _pt_radius <= obj.behind_size[0] + AFFORDACE_PROXIMITY_THRES and
                __pt[1] >= -obj.behind_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.behind_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                _pt_angle >= obj.exist_angle[0]/8 and _pt_angle <= obj.exist_angle[0]):
                return True
            

            return False
        


        elif (isinstance(obj, Middle_Curved_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY', offset_first=True)[0]

            middle_1_mesh_rotation_1 = np.array([0, np.pi, np.pi])
            middle_1_mesh_position_1 = np.array([-obj.middle_size[1], 0, 0])
            middle_1_mesh_rotation_2 = np.array([0, -obj.handle_rotation[1], 0])
            middle_1_mesh_rotation_2_reverse = np.array([0, obj.handle_rotation[1], 0])
            middle_1_mesh_position_1 = adjust_position_from_rotation(middle_1_mesh_position_1, middle_1_mesh_rotation_2)
            middle_1_mesh_position_2 = np.array([
                obj.handle_separation[0] / 2 - obj.front_size[0] / 2 * np.cos(obj.handle_rotation[0]) + obj.front_size[2] * np.sin(obj.handle_rotation[0]) + obj.front_middle_offset[0], 
                obj.front_middle_offset[1],
                -obj.front_size[2] * np.cos(obj.handle_rotation[0]) + obj.front_size[0] / 2 * np.sin(obj.handle_rotation[0])
            ])
            middle_1_mesh_position = middle_1_mesh_position_1 + middle_1_mesh_position_2
            middle_1_mesh_rotation = middle_1_mesh_rotation_1 + middle_1_mesh_rotation_2_reverse

            __pt = apply_transformation(np.array([_pt]), -middle_1_mesh_position, -middle_1_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
            _pt_radius = np.linalg.norm(__pt[[0,2]])
            _pt_angle = np.arctan2(__pt[2], __pt[0])
            
            if (_pt_radius >= obj.middle_size[1] - AFFORDACE_PROXIMITY_THRES and 
                _pt_radius <= obj.middle_size[0] + AFFORDACE_PROXIMITY_THRES and
                __pt[1] >= -obj.middle_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.middle_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                _pt_angle >= obj.exist_angle[0]/10 and _pt_angle <= obj.exist_angle[0]*9/10):
                return True
        

            middle_2_mesh_rotation_1 = np.array([0, np.pi, 0])
            middle_2_mesh_position_1 = np.array([obj.middle_size[1], 0, 0])
            middle_2_mesh_rotation_2 = np.array([0, obj.handle_rotation[1], 0])
            middle_2_mesh_position_1 = adjust_position_from_rotation(middle_2_mesh_position_1, middle_2_mesh_rotation_2)
            middle_2_mesh_position_2 = np.array([
                -obj.handle_separation[0] / 2 + obj.front_size[0] / 2 * np.cos(obj.handle_rotation[0]) - obj.front_size[2] * np.sin(obj.handle_rotation[0]) - obj.front_middle_offset[0], 
                obj.front_middle_offset[1],
                -obj.front_size[2] * np.cos(obj.handle_rotation[0]) + obj.front_size[0] / 2 * np.sin(obj.handle_rotation[0])
            ])
            middle_2_mesh_position = middle_2_mesh_position_1 + middle_2_mesh_position_2
            middle_2_mesh_rotation = middle_2_mesh_rotation_1 + middle_2_mesh_rotation_2

            __pt = apply_transformation(np.array([_pt]), -middle_2_mesh_position, -middle_2_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
            _pt_radius = np.linalg.norm(__pt[[0,2]])
            _pt_angle = np.arctan2(__pt[2], __pt[0])
            
            if (_pt_radius >= obj.middle_size[1] - AFFORDACE_PROXIMITY_THRES and 
                _pt_radius <= obj.middle_size[0] + AFFORDACE_PROXIMITY_THRES and
                __pt[1] >= -obj.middle_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.middle_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                _pt_angle >= obj.exist_angle[0]/10 and _pt_angle <= obj.exist_angle[0]*9/10):
                return True
            

            return False


        
        elif (isinstance(obj, Asymmetric_Straight_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            behind_1_mesh_position_1 = np.array([obj.left_behind_size[0] / 2, 0, -obj.left_behind_size[2] / 2])
            behind_1_mesh_rotation = np.array([0, obj.left_handle_rotation[1], 0])
            behind_1_mesh_position_1 = np.array(adjust_position_from_rotation(behind_1_mesh_position_1, behind_1_mesh_rotation))
            behind_1_mesh_position_2 = np.array([
                -obj.handle_separation[0] / 2 - obj.left_front_size[0] / 2 * np.cos(obj.left_handle_rotation[0]) - obj.left_front_size[2] * np.sin(obj.left_handle_rotation[0]) - obj.left_front_behind_offset[0], 
                obj.left_front_behind_offset[1],
                -obj.left_front_size[2] * np.cos(obj.left_handle_rotation[0]) + obj.left_front_size[0] / 2 * np.sin(obj.left_handle_rotation[0])
            ])
            behind_1_mesh_position = behind_1_mesh_position_1 + behind_1_mesh_position_2

            __pt = inverse_transformation(_pt, behind_1_mesh_position, behind_1_mesh_rotation)

            if (__pt[0] >= -obj.left_behind_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= obj.left_behind_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.left_behind_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.left_behind_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.left_behind_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.left_behind_size[2]/3 + AFFORDACE_PROXIMITY_THRES):
                return True
            

            behind_2_mesh_position_1 = np.array([-obj.right_behind_size[0] / 2, 0, -obj.right_behind_size[2] / 2])
            behind_2_mesh_rotation = np.array([0, -obj.right_handle_rotation[1], 0])
            behind_2_mesh_position_1 = np.array(adjust_position_from_rotation(behind_2_mesh_position_1, behind_2_mesh_rotation))
            behind_2_mesh_position_2 = np.array([
                obj.handle_separation[0] / 2 + obj.right_front_size[0] / 2 * np.cos(obj.right_handle_rotation[0]) + obj.right_front_size[2] * np.sin(obj.right_handle_rotation[0]) + obj.right_front_behind_offset[0], 
                obj.right_front_behind_offset[1] + obj.left_right_offset[0],
                -obj.right_front_size[2] * np.cos(obj.right_handle_rotation[0]) + obj.right_front_size[0] / 2 * np.sin(obj.right_handle_rotation[0])
            ])
            behind_2_mesh_position = behind_2_mesh_position_1 + behind_2_mesh_position_2

            __pt = inverse_transformation(_pt, behind_2_mesh_position, behind_2_mesh_rotation)

            if (__pt[0] >= -obj.right_behind_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= obj.right_behind_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.right_behind_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.right_behind_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.right_behind_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.right_behind_size[2]/3 + AFFORDACE_PROXIMITY_THRES):
                return True
            
            return False


        else:
            return False
    
    
    return is_affordance(obj, pt)






def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
