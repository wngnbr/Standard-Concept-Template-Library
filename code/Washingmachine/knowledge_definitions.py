from concept_template import *
from geometry_template import *
from knowledge_utils import *


def door_affordance(obj, pt):

    def is_affordance(obj, pt):
        if (isinstance(obj, Roller_Door)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            circle_mesh_position = [obj.circle_size[0], 0, obj.circle_size[2] / 2]
            circle_mesh_rotation = [np.pi / 2, 0, 0]
            _pt = inverse_transformation(_pt, circle_mesh_position, circle_mesh_rotation)

            _pt_radius = np.linalg.norm(_pt[[0, 2]])
            _pt_angle = np.arctan2(_pt[2], _pt[0])
            
            return (_pt_radius >= obj.circle_size[1] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.circle_size[0] + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.circle_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.circle_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_angle >= -np.pi/6 and _pt_angle <= np.pi/6)
        

        elif (isinstance(obj, Cuboidal_Door)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = [0, obj.size[1] / 2, 0]
            mesh_rotation = [0, 0, 0]
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)
            
            if (obj.rotation[2] > 0.01):
                return (_pt[0] >= obj.size[0]*2/5 - AFFORDACE_PROXIMITY_THRES and 
                        _pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        _pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                        _pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                        _pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                        _pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES)
            
            elif (obj.rotation[2] < -0.01):
                return (_pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                        _pt[0] <= -obj.size[0]*2/5 + AFFORDACE_PROXIMITY_THRES and 
                        _pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                        _pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                        _pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                        _pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES)
            
            else:
                return (_pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                        _pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        _pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                        _pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                        _pt[2] >= obj.size[2]*2/5 - AFFORDACE_PROXIMITY_THRES and 
                        _pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES)
        

        return False
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
