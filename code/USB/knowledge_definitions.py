from concept_template import *
from geometry_template import *
from knowledge_utils import *


def cap_affordance(obj, pt):

    def is_affordance(obj, pt):
        if (isinstance(obj, Regular_Cap)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = [0, 0, obj.inner_size[2] / 2]
            mesh_rotation = [np.pi / 2, 0, 0]
            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            if (__pt[1] >= -obj.inner_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.inner_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                ((
                    (__pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= -obj.inner_size[0]/2 + obj.inner_outer_offset[0] + AFFORDACE_PROXIMITY_THRES) 
                    or
                    (__pt[0] >= obj.inner_size[0]/2 + obj.inner_outer_offset[0] - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES) 
                ) or 
                (
                    (__pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= -obj.inner_size[1]/2 + obj.inner_outer_offset[1] + AFFORDACE_PROXIMITY_THRES) 
                    or
                    (__pt[2] >= obj.inner_size[1]/2 + obj.inner_outer_offset[1] - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES) 
                ) 
                )):
                return True
            

            mesh_position = [0, 0, (obj.inner_size[2] + obj.size[2]) / 2]
            mesh_rotation = [0, 0, 0]
            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            if (__pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -(obj.size[2] - obj.inner_size[2])/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= (obj.size[2] - obj.inner_size[2])/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            
            return False
        

        elif (isinstance(obj, SquareEnded_Cap)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            radius = obj.proximal_interval[0] / 2 / np.cos(obj.inclination[0]) + obj.size[2] * np.tan(obj.inclination[0])
            back_position = [
                0,
                0,
                -radius * np.sin(obj.inclination[0]) - obj.size[2] * np.cos(obj.inclination[0]) - obj.shaft_offset[0],
            ]
            c_mesh_position = np.array([
                back_position[0] * np.cos(obj.cap_rotation[0]) + back_position[2] * np.sin(obj.cap_rotation[0]),
                back_position[1],
                back_position[0] * np.sin(obj.cap_rotation[0]) + back_position[2] * np.cos(obj.cap_rotation[0]) + obj.shaft_offset[0],
            ])
            c_mesh_rotation = np.array([0, np.pi, np.pi / 2])
            _pt = apply_transformation(np.array([_pt]), -np.array(c_mesh_position), -np.array(c_mesh_rotation), rotation_order='ZYX', offset_first=True)[0]

            _pt_radius = np.linalg.norm(_pt[[0, 2]])
            _pt_angle = np.arctan2(_pt[2], _pt[0])
            
            return (_pt_radius >= radius - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= radius + obj.size[1] + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES)
                    # _pt_angle >= 0 and _pt_angle <= np.pi + obj.inclination[0] * 2)
        

        elif (isinstance(obj, RoundEnded_Cap)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            radius = obj.proximal_interval[0] / 2 / np.cos(obj.inclination[0]) + obj.size[2] * np.tan(obj.inclination[0])
            back_position = [
                0,
                0,
                -radius * np.sin(obj.inclination[0]) - obj.size[2] * np.cos(obj.inclination[0]) - obj.shaft_offset[0],
            ]
            c_mesh_position = np.array([
                back_position[0] * np.cos(obj.cap_rotation[0]) + back_position[2] * np.sin(obj.cap_rotation[0]),
                back_position[1],
                back_position[0] * np.sin(obj.cap_rotation[0]) + back_position[2] * np.cos(obj.cap_rotation[0]) + obj.shaft_offset[0],
            ])
            c_mesh_rotation = np.array([0, np.pi, np.pi / 2])
            _pt = apply_transformation(np.array([_pt]), -np.array(c_mesh_position), -np.array(c_mesh_rotation), rotation_order='ZYX', offset_first=True)[0]

            _pt_radius = np.linalg.norm(_pt[[0, 2]])
            _pt_angle = np.arctan2(_pt[2], _pt[0])
            
            return (_pt_radius >= radius - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= radius + obj.size[1] + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES)
                    # _pt_angle >= 0 and _pt_angle <= np.pi + obj.inclination[0] * 2)
        

        return False
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
