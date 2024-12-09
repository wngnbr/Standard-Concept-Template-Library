from concept_template import *
from geometry_template import *
from knowledge_utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Cuboidal_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ', offset_first=True)[0]

            mesh_position = np.array([0, 0, obj.size[2] / 2])
            mesh_rotation = np.array([0, 0, 0])

            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            return (__pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Trifold_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ', offset_first=True)[0]
            
            mesh_position = np.array([0, 0, obj.mounting_size[2] + obj.grip_size[2] / 2])
            mesh_rotation = np.array([0, 0, 0])
            
            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            return (__pt[0] >= -obj.grip_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[0] <= obj.grip_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.grip_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.grip_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] >= -obj.grip_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] <= obj.grip_size[2]/2 + AFFORDACE_PROXIMITY_THRES )
    


        elif (isinstance(obj, Trifold_Curve_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ', offset_first=True)[0]

            curve_z_offset = obj.mounting_size[2] - np.sqrt(obj.curve_size[1]**2 - (obj.mounting_seperation[0] / 2)**2)
            mesh_position = np.array([0, 0, curve_z_offset])
            mesh_rotation = np.array([0, -np.pi / 2 + obj.curve_exist_angle[0] / 2, np.pi / 2])

            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            
            _pt_radius = np.linalg.norm(__pt[[0,2]])
            _pt_angle = np.arccos(__pt[0] / np.sqrt(__pt[0]**2 + __pt[2]**2))
            
            return (_pt_radius >= obj.curve_size[1] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.curve_size[0] + AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.curve_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.curve_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt_angle >= obj.curve_exist_angle[0]/3 and _pt_angle <= obj.curve_exist_angle[0]*2/3)
                

        elif (isinstance(obj, Curve_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ', offset_first=True)[0]

            curve_z_offset = -obj.curve_size[0] * np.cos(obj.curve_exist_angle[0] / 2)
            mesh_position = np.array([0, 0, curve_z_offset])
            mesh_rotation = np.array([0, -np.pi / 2 + obj.curve_exist_angle[0] / 2, np.pi / 2])

            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            
            _pt_radius = np.linalg.norm(__pt[[0,2]])
            _pt_angle = np.arccos(__pt[0] / np.sqrt(__pt[0]**2 + __pt[2]**2))
            
            return (_pt_radius >= obj.curve_size[1] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.curve_size[0] + AFFORDACE_PROXIMITY_THRES and
                    abs(__pt[1]) <= obj.curve_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_angle >= obj.curve_exist_angle[0]/3 and _pt_angle <= obj.curve_exist_angle[0]*2/3)
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
