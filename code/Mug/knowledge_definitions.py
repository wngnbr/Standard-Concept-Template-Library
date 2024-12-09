from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Trifold_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY', offset_first=True)[0]

            vertical_y_offset = (-obj.horizontal_length[0] * np.sin(obj.horizontal_rotation[0]) - obj.horizontal_length[1] * np.sin(obj.horizontal_rotation[1])) / 2
            vertical_z_offset = (obj.horizontal_length[1] * np.cos(obj.horizontal_rotation[1]) + obj.mounting_offset[0] + obj.horizontal_length[0] * np.cos(obj.horizontal_rotation[0])) / 2
            delta_y = obj.horizontal_separation[0] - obj.horizontal_length[0] * np.sin(obj.horizontal_rotation[0]) + obj.horizontal_length[1] * np.sin(obj.horizontal_rotation[1])
            delta_z = obj.mounting_offset[0] - obj.horizontal_length[1] * np.cos(obj.horizontal_rotation[1]) + obj.horizontal_length[0] * np.cos(obj.horizontal_rotation[0])
            vertical_rotation = np.arctan(delta_z / delta_y)
            vertical_length = np.sqrt(delta_y * delta_y + delta_z * delta_z) + obj.horizontal_thickness[1]
            vertical_mesh_position = np.array([
                0, 
                vertical_y_offset, 
                vertical_z_offset + obj.vertical_thickness[1] / 2
            ])
            vertical_mesh_rotation = np.array([-vertical_rotation, 0, 0])
            _pt = inverse_transformation(_pt, vertical_mesh_position, vertical_mesh_rotation)

            return (_pt[0] >= -obj.vertical_thickness[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[0] <= obj.vertical_thickness[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -vertical_length/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[1] <= vertical_length/2 + AFFORDACE_PROXIMITY_THRES and
                    _pt[2] >= -obj.vertical_thickness[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[2] <= obj.vertical_thickness[1]/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Curved_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = np.array([0, 0, 0])
            mesh_rotation = np.array([0, 0, -np.pi / 2])
            _pt = apply_transformation(np.array([_pt]), -mesh_position, -mesh_rotation, rotation_order='ZYX', offset_first=True)[0]

            _pt_radius = np.linalg.norm(_pt[[0, 2]])
            _pt_angle = np.arctan2(_pt[2], _pt[0])

            return (_pt_radius <= obj.radius[0] + obj.radius[1] + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius >= obj.radius[0] - obj.radius[1] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_angle >= obj.central_angle[0]/3 and _pt_angle <= obj.central_angle[0]*2/3)
        
        else:
            return False
    
    
    return is_affordance(obj, pt)



def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
