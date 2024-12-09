from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Cuboidal_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY', offset_first=True)[0]

            return (_pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[1] <= obj.size[1]/3 + AFFORDACE_PROXIMITY_THRES and
                    _pt[2] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[2] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES)


        elif (isinstance(obj, T_Shaped_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY', offset_first=True)[0]


            return (_pt[0] >= -obj.main_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[0] <= obj.main_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.main_size[1]/3 - AFFORDACE_PROXIMITY_THRES and
                    _pt[1] <= obj.main_size[1]/3 + AFFORDACE_PROXIMITY_THRES and
                    _pt[2] >= -obj.main_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[2] <= obj.main_size[1]/2 + AFFORDACE_PROXIMITY_THRES)
        


        elif (isinstance(obj, Cylindrical_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            return (_pt[1] >= -obj.size[2]/3 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.size[2]/3 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= (obj.size[0] + obj.size[1]) / 2 - _pt[1] * (obj.size[1] - obj.size[0]) / obj.size[2] + AFFORDACE_PROXIMITY_THRES)


        elif (isinstance(obj, Curved_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY', offset_first=True)[0]

            mesh_position = np.array([0, 0, (obj.radius[0] + obj.radius[1]) / 2])
            mesh_rotation = np.array([np.pi / 2, np.pi / 2, 0])
            _pt = apply_transformation(np.array([_pt]), -mesh_position, -mesh_rotation, rotation_order='ZYX', offset_first=True)[0]

            _pt_radius = np.linalg.norm(_pt[[0, 2]])
            _pt_angle = np.arctan2(_pt[2], _pt[0])

            return (_pt[1] >= -obj.thickness[0]/3 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.thickness[0]/3 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.radius[0] + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius >= obj.radius[1] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_angle >= obj.exist_angle[0] / 9 and 
                    _pt_angle <= obj.exist_angle[0] * 8 / 9)
        


        elif (isinstance(obj, Multideck_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            left_mesh_position = np.array([
                obj.beside_seperation[0] / 2 + obj.beside_offset[0],
                obj.bottom_size[1] / 2,
                obj.beside_offset[1]
            ])
            left_mesh_rotation = np.array([0, 0, 0])
            __pt = apply_transformation(np.array([_pt]), -left_mesh_position, -left_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]

            if (__pt[0] >= -obj.beside_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[0] <= obj.beside_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.beside_size[1]/3 - AFFORDACE_PROXIMITY_THRES and
                __pt[1] <= obj.beside_size[1]/3 + AFFORDACE_PROXIMITY_THRES and
                __pt[2] >= -obj.beside_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[2] <= obj.beside_size[1]/2 + AFFORDACE_PROXIMITY_THRES):
                return True

            right_mesh_position = np.array([
                -obj.beside_seperation[0] / 2 + obj.beside_offset[0],
                obj.bottom_size[1] / 2,
                obj.beside_offset[1]
            ])
            right_mesh_rotation = np.array([0, 0, 0])
            __pt = apply_transformation(np.array([_pt]), -right_mesh_position, -right_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]

            if (__pt[0] >= -obj.beside_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[0] <= obj.beside_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.beside_size[1]/3 - AFFORDACE_PROXIMITY_THRES and
                __pt[1] <= obj.beside_size[1]/3 + AFFORDACE_PROXIMITY_THRES and
                __pt[2] >= -obj.beside_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[2] <= obj.beside_size[1]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            
            return False
        


        elif (isinstance(obj, Enveloping_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            left_mesh_position = np.array([
                obj.thickness[1] / 2,
                obj.thickness[0] / 2,
                (obj.size[2] - obj.thickness[2]) / 2
            ])
            left_mesh_rotation = np.array([0, 0, 0])
            __pt = apply_transformation(np.array([_pt]), -left_mesh_position, -left_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]

            if (__pt[0] >= -(obj.size[0] - obj.thickness[1])/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[0] <= (obj.size[0] - obj.thickness[1])/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -(obj.size[1] - obj.thickness[0])/3 - AFFORDACE_PROXIMITY_THRES and
                __pt[1] <= (obj.size[1] - obj.thickness[0])/3 + AFFORDACE_PROXIMITY_THRES and
                __pt[2] >= -obj.thickness[2]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[2] <= obj.thickness[2]/2 + AFFORDACE_PROXIMITY_THRES):
                return True

            right_mesh_position = np.array([
                obj.thickness[1] / 2,
                obj.thickness[0] / 2,
                -(obj.size[2] - obj.thickness[2]) / 2
            ])
            right_mesh_rotation = np.array([0, 0, 0])
            __pt = apply_transformation(np.array([_pt]), -right_mesh_position, -right_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]

            if (__pt[0] >= -(obj.size[0] - obj.thickness[1])/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[0] <= (obj.size[0] - obj.thickness[1])/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -(obj.size[1] - obj.thickness[0])/3 - AFFORDACE_PROXIMITY_THRES and
                __pt[1] <= (obj.size[1] - obj.thickness[0])/3 + AFFORDACE_PROXIMITY_THRES and
                __pt[2] >= -obj.thickness[2]/2 - AFFORDACE_PROXIMITY_THRES and
                __pt[2] <= obj.thickness[2]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            
            return False
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
