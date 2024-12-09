from concept_template import *
from geometry_template import *
from knowledge_utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Ring_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            arm_mesh_position_1 = np.array([obj.arm_radius[0], 0, 0])
            arm_mesh_rotation = np.array([0, -obj.arm_rotation[0], 0])
            arm_mesh_position_1 = adjust_position_from_rotation(arm_mesh_position_1, arm_mesh_rotation)
            arm_mesh_position_2 = np.array([
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ])
            arm_mesh_position = arm_mesh_position_1 + arm_mesh_position_2

            __pt = inverse_transformation(_pt, arm_mesh_position, arm_mesh_rotation)
            _pt_radius = np.linalg.norm(__pt[[0, 2]])

            return (_pt_radius >= obj.arm_radius[3] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.arm_radius[1] + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.arm_radius[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.arm_radius[0]/2 - AFFORDACE_PROXIMITY_THRES)

    
                

        elif (isinstance(obj, Half_Ring_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            arm_mesh_position_1 = np.array([obj.arm_radius[0], 0, 0])
            arm_mesh_rotation = np.array([0, -obj.arm_rotation[0], 0])
            arm_mesh_position_1 = adjust_position_from_rotation(arm_mesh_position_1, arm_mesh_rotation)
            arm_mesh_position_2 = np.array([
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ])
            arm_mesh_position = arm_mesh_position_1 + arm_mesh_position_2

            __pt = inverse_transformation(_pt, arm_mesh_position, arm_mesh_rotation)
            _pt_radius = np.linalg.norm(__pt[[0, 2]])
            _pt_angle = np.arctan2(__pt[2], __pt[0])

            if (_pt_radius >= obj.arm_radius[3] - AFFORDACE_PROXIMITY_THRES and 
                _pt_radius <= obj.arm_radius[1] + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and
                _pt_angle >= 0 and _pt_angle <= np.pi):
                return True


            horizontal_mesh_position_1 = np.array([
                obj.arm_radius[0], 
                0,
                -obj.arm_horizontal_thickness[0] / 2
            ])
            horizontal_mesh_rotation = np.array([0, -obj.arm_rotation[0], 0])
            horizontal_mesh_position_1 = adjust_position_from_rotation(horizontal_mesh_position_1, horizontal_mesh_rotation)
            horizontal_mesh_position_2 = np.array([
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ])
            horizontal_mesh_position = horizontal_mesh_position_1 + horizontal_mesh_position_2
            __pt = inverse_transformation(_pt, horizontal_mesh_position, horizontal_mesh_rotation)

            if (__pt[0] >= -obj.arm_radius[0] - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= obj.arm_radius[0] + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.arm_horizontal_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.arm_horizontal_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            
            return False
        


        elif (isinstance(obj, Double_Curved_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            bottom_mesh_position_1 = np.array([obj.arm_radius[1], 0, 0])
            bottom_mesh_rotation = np.array([0, -obj.arm_rotation[0], 0])
            bottom_mesh_position_1 = adjust_position_from_rotation(bottom_mesh_position_1, bottom_mesh_rotation)
            bottom_mesh_position_2 = [
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ]
            bottom_mesh_position = list_add(bottom_mesh_position_1, bottom_mesh_position_2)
            bottom_mesh_rotation[1] += -np.pi + obj.arm_exist_angle[1] / 2

            __pt = inverse_transformation(_pt, bottom_mesh_position, bottom_mesh_rotation)
            _pt_radius = np.linalg.norm(__pt[[0, 2]])
            _pt_angle = np.arctan2(__pt[2], __pt[0])

            if (_pt_radius >= obj.arm_radius[1] - obj.arm_thickness[1] - AFFORDACE_PROXIMITY_THRES and 
                _pt_radius <= obj.arm_radius[1] + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and
                _pt_angle >= 0 and _pt_angle <= obj.arm_exist_angle[1]):
                return True
            

            top_mesh_position_1 = np.array([
                -obj.arm_radius[0] * np.cos(obj.arm_exist_angle[0] / 2) + obj.arm_radius[1] * (1 - np.cos(obj.arm_exist_angle[1] / 2)) + obj.arm_top_bottom_separation[0], 
                0, 
                0
            ])
            top_mesh_rotation = np.array([0, -obj.arm_rotation[0], 0])
            top_mesh_position_1 = adjust_position_from_rotation(top_mesh_position_1, top_mesh_rotation)
            top_mesh_position_2 = [
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ]
            top_mesh_position = list_add(top_mesh_position_1, top_mesh_position_2)
            top_mesh_rotation[1] += obj.arm_exist_angle[1] / 2

            __pt = inverse_transformation(_pt, top_mesh_position, top_mesh_rotation)
            _pt_radius = np.linalg.norm(__pt[[0, 2]])
            _pt_angle = np.arctan2(__pt[2], __pt[0])

            if (_pt_radius >= obj.arm_radius[0] - obj.arm_thickness[1] - AFFORDACE_PROXIMITY_THRES and 
                _pt_radius <= obj.arm_radius[0] + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and
                _pt_angle >= 0 and _pt_angle <= obj.arm_exist_angle[0]):
                return True
            

            arm_delta_z = obj.arm_radius[0] * np.sin(obj.arm_exist_angle[0] / 2) - obj.arm_radius[1] * np.sin(obj.arm_exist_angle[1] / 2)
            box_length = np.sqrt(arm_delta_z * arm_delta_z + obj.arm_top_bottom_separation[0] * obj.arm_top_bottom_separation[0])
            box_rotation = np.arctan(arm_delta_z / obj.arm_top_bottom_separation[0])
            left_mesh_position_1 = np.array([
                obj.arm_radius[1] * (1 - np.cos(obj.arm_exist_angle[1] / 2)) + obj.arm_top_bottom_separation[0] / 2 + obj.arm_thickness[1] / 2 * np.sin(box_rotation), 
                0,
                -(obj.arm_radius[0] * np.sin(obj.arm_exist_angle[0] / 2) + obj.arm_radius[1] * np.sin(obj.arm_exist_angle[1] / 2)) / 2 + obj.arm_thickness[1] / 2 * np.cos(box_rotation)
            ])
            left_mesh_rotation = np.array([0, -obj.arm_rotation[0], 0])
            left_mesh_position_1 = adjust_position_from_rotation(left_mesh_position_1, left_mesh_rotation)

            left_mesh_position_2 = np.array([
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ])
            left_mesh_position = left_mesh_position_1 + left_mesh_position_2
            left_mesh_rotation[1] += box_rotation

            __pt = inverse_transformation(_pt, left_mesh_position, left_mesh_rotation)

            if (__pt[0] >= -box_length/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= box_length/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.arm_thickness[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.arm_thickness[1]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            

            right_mesh_position_1 = np.array([
                obj.arm_radius[1] * (1 - np.cos(obj.arm_exist_angle[1] / 2)) + obj.arm_top_bottom_separation[0] / 2 + obj.arm_thickness[1] / 2 * np.sin(box_rotation), 
                0,
                (obj.arm_radius[0] * np.sin(obj.arm_exist_angle[0] / 2) + obj.arm_radius[1] * np.sin(obj.arm_exist_angle[1] / 2)) / 2 - obj.arm_thickness[1] / 2 * np.cos(box_rotation)
            ])
            right_mesh_rotation = np.array([0, -obj.arm_rotation[0], 0])
            right_mesh_position_1 = adjust_position_from_rotation(right_mesh_position_1, right_mesh_rotation)

            right_mesh_position_2 = np.array([
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ])
            right_mesh_position = right_mesh_position_1 + right_mesh_position_2
            right_mesh_rotation[1] += -box_rotation

            __pt = inverse_transformation(_pt, right_mesh_position, right_mesh_rotation)

            if (__pt[0] >= -box_length/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= box_length/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.arm_thickness[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.arm_thickness[1]/2 + AFFORDACE_PROXIMITY_THRES):
                return True

            return False


        
        elif (isinstance(obj, Triple_Curved_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            bottom_mesh_position_1 = np.array([obj.arm_radius[0], 0, 0])
            bottom_mesh_rotation = np.array([0, -obj.arm_rotation[0], 0])
            bottom_mesh_position_1 = adjust_position_from_rotation(bottom_mesh_position_1, bottom_mesh_rotation)
            bottom_mesh_position_2 = [
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ]
            bottom_mesh_position = list_add(bottom_mesh_position_1, bottom_mesh_position_2)
            bottom_mesh_rotation[1] += -np.pi + obj.arm_exist_angle[0] / 2

            __pt = inverse_transformation(_pt, bottom_mesh_position, bottom_mesh_rotation)
            _pt_radius = np.linalg.norm(__pt[[0, 2]])
            _pt_angle = np.arctan2(__pt[2], __pt[0])

            if (_pt_radius >= obj.arm_radius[0] - obj.arm_thickness[1] - AFFORDACE_PROXIMITY_THRES and 
                _pt_radius <= obj.arm_radius[0] + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and
                _pt_angle >= 0 and _pt_angle <= obj.arm_exist_angle[0]):
                return True
            

            top_mesh_position_1 = np.array([
                -obj.arm_radius[0] * np.cos(obj.arm_exist_angle[0] / 2) + obj.arm_radius[0] * (1 - np.cos(obj.arm_exist_angle[0] / 2)) + obj.arm_top_bottom_separation[0], 
                0, 
                0
            ])
            top_mesh_rotation = np.array([0, -obj.arm_rotation[0], 0])
            top_mesh_position_1 = adjust_position_from_rotation(top_mesh_position_1, top_mesh_rotation)
            top_mesh_position_2 = [
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ]
            top_mesh_position = list_add(top_mesh_position_1, top_mesh_position_2)
            top_mesh_rotation[1] += obj.arm_exist_angle[0] / 2

            __pt = inverse_transformation(_pt, top_mesh_position, top_mesh_rotation)
            _pt_radius = np.linalg.norm(__pt[[0, 2]])
            _pt_angle = np.arctan2(__pt[2], __pt[0])

            if (_pt_radius >= obj.arm_radius[0] - obj.arm_thickness[1] - AFFORDACE_PROXIMITY_THRES and 
                _pt_radius <= obj.arm_radius[0] + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and
                _pt_angle >= 0 and _pt_angle <= obj.arm_exist_angle[0]):
                return True
            

            arm_delta_z = obj.arm_radius[0] * np.sin(obj.arm_exist_angle[0] / 2) - obj.arm_radius[0] * np.sin(obj.arm_exist_angle[0] / 2)
            box_length = np.sqrt(arm_delta_z * arm_delta_z + obj.arm_top_bottom_separation[0] * obj.arm_top_bottom_separation[0])
            box_rotation = np.arctan(arm_delta_z / obj.arm_top_bottom_separation[0])
            left_mesh_position_1 = np.array([
                obj.arm_radius[0] * (1 - np.cos(obj.arm_exist_angle[0] / 2)) + obj.arm_top_bottom_separation[0] / 2 + obj.arm_thickness[1] / 2 * np.sin(box_rotation), 
                0,
                -(obj.arm_radius[0] * np.sin(obj.arm_exist_angle[0] / 2) + obj.arm_radius[0] * np.sin(obj.arm_exist_angle[0] / 2)) / 2 + obj.arm_thickness[1] / 2 * np.cos(box_rotation)
            ])
            left_mesh_rotation = np.array([0, -obj.arm_rotation[0], 0])
            left_mesh_position_1 = adjust_position_from_rotation(left_mesh_position_1, left_mesh_rotation)

            left_mesh_position_2 = np.array([
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ])
            left_mesh_position = left_mesh_position_1 + left_mesh_position_2
            left_mesh_rotation[1] += box_rotation

            __pt = inverse_transformation(_pt, left_mesh_position, left_mesh_rotation)

            if (__pt[0] >= -box_length/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= box_length/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.arm_thickness[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.arm_thickness[1]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            

            right_mesh_position_1 = np.array([
                obj.arm_radius[0] * (1 - np.cos(obj.arm_exist_angle[0] / 2)) + obj.arm_top_bottom_separation[0] / 2 + obj.arm_thickness[1] / 2 * np.sin(box_rotation), 
                0,
                (obj.arm_radius[0] * np.sin(obj.arm_exist_angle[0] / 2) + obj.arm_radius[0] * np.sin(obj.arm_exist_angle[0] / 2)) / 2 - obj.arm_thickness[1] / 2 * np.cos(box_rotation)
            ])
            right_mesh_rotation = np.array([0, -obj.arm_rotation[0], 0])
            right_mesh_position_1 = adjust_position_from_rotation(right_mesh_position_1, right_mesh_rotation)

            right_mesh_position_2 = np.array([
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ])
            right_mesh_position = right_mesh_position_1 + right_mesh_position_2
            right_mesh_rotation[1] += -box_rotation

            __pt = inverse_transformation(_pt, right_mesh_position, right_mesh_rotation)

            if (__pt[0] >= -box_length/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= box_length/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.arm_thickness[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.arm_thickness[1]/2 + AFFORDACE_PROXIMITY_THRES):
                return True

            return False
        

        elif (isinstance(obj, Cuboidal_Ring_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            bottom_mesh_position_1 = np.array([0, 0, 0])
            bottom_mesh_rotation = np.array([0, -obj.arm_rotation[0], 0])
            bottom_mesh_position_1 = adjust_position_from_rotation(bottom_mesh_position_1, bottom_mesh_rotation)
            bottom_mesh_position_2 = [
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ]
            bottom_mesh_position = list_add(bottom_mesh_position_1, bottom_mesh_position_2)
            bottom_mesh_rotation[1] += obj.arm_top_bottom_rotation[1]

            __pt = inverse_transformation(_pt, bottom_mesh_position, bottom_mesh_rotation)
            _pt_radius = np.linalg.norm(__pt[[0, 2]])
            _pt_angle = np.arctan2(__pt[2], __pt[0])

            if (__pt[0] >= -obj.arm_thickness[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= obj.arm_thickness[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.arm_length[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.arm_length[1]/2 + AFFORDACE_PROXIMITY_THRES):
                return True


            top_mesh_position_1 = np.array([
                obj.arm_top_bottom_separation[0], 
                0, 
                obj.arm_top_bottom_offset[0]
            ])
            top_mesh_rotation = np.array([0, -obj.arm_rotation[0], 0])
            top_mesh_position_1 = adjust_position_from_rotation(top_mesh_position_1, top_mesh_rotation)
            top_mesh_position_2 = [
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ]
            top_mesh_position = list_add(top_mesh_position_1, top_mesh_position_2)
            top_mesh_rotation[1] += obj.arm_top_bottom_rotation[0]

            __pt = inverse_transformation(_pt, top_mesh_position, top_mesh_rotation)
            _pt_radius = np.linalg.norm(__pt[[0, 2]])
            _pt_angle = np.arctan2(__pt[2], __pt[0])

            if (__pt[0] >= -obj.arm_thickness[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= obj.arm_thickness[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.arm_length[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.arm_length[0]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            

            left_arm_length_x = obj.arm_top_bottom_separation[0] - obj.arm_length[0] * np.sin(obj.arm_top_bottom_rotation[0]) / 2 + obj.arm_length[1] * np.sin(obj.arm_top_bottom_rotation[1]) / 2
            left_arm_length_z = obj.arm_length[0] * np.cos(obj.arm_top_bottom_rotation[0]) / 2 - obj.arm_length[1] * np.cos(obj.arm_top_bottom_rotation[1]) / 2 - obj.arm_top_bottom_offset[0]
            left_arm_length_t = np.sqrt(left_arm_length_x * left_arm_length_x + left_arm_length_z * left_arm_length_z)
            left_arm_offset_x = (obj.arm_top_bottom_separation[0] - obj.arm_length[0] * np.sin(obj.arm_top_bottom_rotation[0]) / 2 - obj.arm_length[1] * np.sin(obj.arm_top_bottom_rotation[1]) / 2) / 2
            left_arm_offset_z = (obj.arm_top_bottom_offset[0] - obj.arm_length[0] * np.cos(obj.arm_top_bottom_rotation[0]) / 2 - obj.arm_length[1] * np.cos(obj.arm_top_bottom_rotation[1]) / 2) / 2
            left_arm_rotation = np.arctan(left_arm_length_z / left_arm_length_x)
            left_mesh_position_1 = [
                left_arm_offset_x, 
                0,
                left_arm_offset_z
            ]
            left_mesh_rotation = [0, -obj.arm_rotation[0], 0]
            left_mesh_position_1 = adjust_position_from_rotation(left_mesh_position_1, left_mesh_rotation)
            left_mesh_position_2 = [
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ]
            left_mesh_position = list_add(left_mesh_position_1, left_mesh_position_2)
            left_mesh_rotation[1] += left_arm_rotation

            __pt = inverse_transformation(_pt, left_mesh_position, left_mesh_rotation)
            if (__pt[0] >= -left_arm_length_t/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= left_arm_length_t/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.arm_length[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.arm_length[1]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            

            right_arm_length_x = obj.arm_top_bottom_separation[0] + obj.arm_length[0] * np.sin(obj.arm_top_bottom_rotation[0]) / 2 - obj.arm_length[1] * np.sin(obj.arm_top_bottom_rotation[1]) / 2
            right_arm_length_z = obj.arm_length[0] * np.cos(obj.arm_top_bottom_rotation[0]) / 2 - obj.arm_length[1] * np.cos(obj.arm_top_bottom_rotation[1]) / 2 + obj.arm_top_bottom_offset[0]
            right_arm_length_t = np.sqrt(right_arm_length_x * right_arm_length_x + right_arm_length_z * right_arm_length_z)
            right_arm_offset_x = (obj.arm_top_bottom_separation[0] + obj.arm_length[0] * np.sin(obj.arm_top_bottom_rotation[0]) / 2 + obj.arm_length[1] * np.sin(obj.arm_top_bottom_rotation[1]) / 2) / 2
            right_arm_offset_z = (obj.arm_top_bottom_offset[0] + obj.arm_length[0] * np.cos(obj.arm_top_bottom_rotation[0]) / 2 + obj.arm_length[1] * np.cos(obj.arm_top_bottom_rotation[1]) / 2) / 2
            right_arm_rotation = np.arctan(right_arm_length_z / right_arm_length_x)
            right_mesh_position_1 = [
                right_arm_offset_x, 
                0,
                right_arm_offset_z
            ]
            right_mesh_rotation = [0, -obj.arm_rotation[0], 0]
            right_mesh_position_1 = adjust_position_from_rotation(right_mesh_position_1, right_mesh_rotation)
            right_mesh_position_2 = [
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ]
            right_mesh_position = list_add(right_mesh_position_1, right_mesh_position_2)
            right_mesh_rotation[1] += -right_arm_rotation

            __pt = inverse_transformation(_pt, right_mesh_position, right_mesh_rotation)
            if (__pt[0] >= -right_arm_length_t/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= right_arm_length_t/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.arm_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.arm_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.arm_length[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.arm_length[1]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            
            return False
        

        elif (isinstance(obj, Cuboidal_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            arm_mesh_position_1 = [obj.arm_size[0] / 2, 0, 0]
            arm_mesh_rotation = [0, -obj.arm_rotation[0], 0]
            arm_mesh_position_1 = adjust_position_from_rotation(arm_mesh_position_1, arm_mesh_rotation)
            arm_mesh_position_2 = [
                obj.root_size[0] * np.cos(obj.root_rotation[0]) - obj.arm_offset[0] * np.tan(obj.root_rotation[0]), 
                0,
                obj.root_size[0] * np.sin(obj.root_rotation[0]) + obj.arm_offset[0]
            ]
            arm_mesh_position = list_add(arm_mesh_position_1, arm_mesh_position_2)
            _pt = inverse_transformation(_pt, arm_mesh_position, arm_mesh_rotation)

            return (_pt[0] >= -obj.arm_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[0] <= obj.arm_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.arm_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.arm_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[2] >= -obj.arm_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[2] <= obj.arm_size[2]/2 + AFFORDACE_PROXIMITY_THRES)


        else:
            return False
    
    
    return is_affordance(obj, pt)




def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
