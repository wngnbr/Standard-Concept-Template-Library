from concept_template import *
from geometry_template import *
from knowledge_utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Trifold_Clip)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = np.array([
                0,
                -obj.clip_root_size[1] / 2,
                obj.clip_root_size[2] / 2
            ])
            mesh_rotation = np.array([0, 0, 0])

            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            if (__pt[0] >= -obj.clip_root_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= obj.clip_root_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.clip_root_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.clip_root_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.clip_root_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.clip_root_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            
            mesh_position = np.array([
                0,
                -obj.clip_vertical_size[1] / 2 - obj.clip_offset[0],
                obj.clip_root_size[2] + obj.clip_vertical_size[2] / 2
            ])
            mesh_rotation = np.array([0, 0, 0])

            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            if (__pt[0] >= -obj.clip_vertical_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= obj.clip_vertical_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.clip_vertical_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.clip_vertical_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.clip_vertical_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.clip_vertical_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            
            mesh_position = np.array([
                0,
                -obj.clip_vertical_size[1] - obj.clip_offset[0] + obj.clip_tip_size[1] / 2 - obj.clip_offset[1],
                obj.clip_root_size[2] - obj.clip_tip_size[2] / 2
            ])
            mesh_rotation = np.array([0, 0, 0])

            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            if (__pt[0] >= -obj.clip_tip_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= obj.clip_tip_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.clip_tip_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.clip_tip_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.clip_tip_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.clip_tip_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            
            return False
    
                

        elif (isinstance(obj, Curved_Clip)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY', offset_first=True)[0]

            mesh_position = np.array([
                0, 
                -obj.clip_curve_size[0] * np.sin(obj.clip_curve_exist_angle[0] / 2), 
                -obj.clip_curve_size[0] * np.cos(obj.clip_curve_exist_angle[0] / 2)
            ])
            mesh_rotation = np.array([
                0,
                -np.pi / 2 + obj.clip_curve_exist_angle[0] / 2,
                np.pi / 2
            ])

            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            
            _pt_radius = np.linalg.norm(__pt[[0,2]])
            _pt_angle = np.arctan2(__pt[2], __pt[0])
            
            return (_pt_radius >= obj.clip_curve_size[1] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.clip_curve_size[0] + AFFORDACE_PROXIMITY_THRES and
                    __pt[1] >= -obj.clip_curve_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.clip_curve_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_angle >= 0 and _pt_angle <= obj.clip_curve_exist_angle[0])
        
        else:
            return False
    
    
    return is_affordance(obj, pt)




def button_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Cylindrical_Button)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = np.array([0, obj.size[2]/2, 0])
            mesh_rotation = np.array([0, 0, 0])
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            return (_pt[1] >= - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= (obj.size[0] + obj.size[1]) / 2 - _pt[1] * (obj.size[1] - obj.size[0]) / obj.size[2] + AFFORDACE_PROXIMITY_THRES)
    
                

        elif (isinstance(obj, Bistratal_Button)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = np.array([0, obj.bottom_size[2] + obj.top_size[2] / 2, 0])
            mesh_rotation = np.array([0, 0, 0])
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            return (_pt[1] >= - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.top_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= (obj.top_size[0] + obj.top_size[1]) / 2 - _pt[1] * (obj.top_size[1] - obj.top_size[0]) / obj.top_size[2] + AFFORDACE_PROXIMITY_THRES)
        
        else:
            return False
    
    
    return is_affordance(obj, pt)




def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
