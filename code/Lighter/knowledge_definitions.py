from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def wheel_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Simplified_Wheel)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = np.array([0, 0, 0])
            mesh_rotation = np.array([0, 0, np.pi / 2])
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            return (_pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.size[0] + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Standard_Wheel)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            left_mesh_position = np.array([
                (obj.middle_size[1] + obj.beside_size[1]) / 2, 
                0,
                0
            ])
            left_mesh_rotation = np.array([0, 0, np.pi / 2])
            __pt = inverse_transformation(_pt, left_mesh_position, left_mesh_rotation)

            _pt_radius = np.linalg.norm(__pt[[0, 2]])

            if (__pt[1] >= -obj.beside_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.beside_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                _pt_radius <= obj.beside_size[0] + AFFORDACE_PROXIMITY_THRES):
                return True
            

            right_mesh_position = np.array([
                -(obj.middle_size[1] + obj.beside_size[1]) / 2, 
                0,
                0
            ])
            right_mesh_rotation = np.array([0, 0, np.pi / 2])
            __pt = inverse_transformation(_pt, right_mesh_position, right_mesh_rotation)

            _pt_radius = np.linalg.norm(__pt[[0, 2]])

            if (__pt[1] >= -obj.beside_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.beside_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                _pt_radius <= obj.beside_size[0] + AFFORDACE_PROXIMITY_THRES):
                return True

            return False
        
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def button_affordance(obj, pt):

    def is_affordance(obj, pt):
        if (isinstance(obj, L_Shaped_Button)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = np.array([
                0, 
                obj.bottom_size[1] + obj.top_size[1] / 2,
                obj.top_bottom_offset[0]
            ])
            mesh_rotation = np.array([0, 0, 0])
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            return (_pt[0] >= -obj.top_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[0] <= obj.top_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.top_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[1] <= obj.top_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    _pt[2] >= -obj.top_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[2] <= obj.top_size[1]/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Double_Cambered_Button)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = np.array([
                0, 
                obj.bottom_size[1] + obj.top_size[1] / 2,
                -obj.top_size[2] / 2 + obj.top_bottom_offset[0]
            ])
            mesh_rotation = np.array([0, np.pi, 0])
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            return (_pt[1] >= -obj.top_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.top_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.top_size[0] / 2 + AFFORDACE_PROXIMITY_THRES)
        
        else:
            return False

    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
