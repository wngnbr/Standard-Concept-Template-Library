from concept_template import *
from geometry_template import *
from knowledge_utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):
        
        if (isinstance(obj, Trifold_Handle)):
            delta_x = obj.vertical_separation[0] - obj.vertical_length[0] * np.sin(obj.vertical_rotation[0]) + obj.vertical_length[1] * np.sin(obj.vertical_rotation[1])
            delta_y = obj.vertical_length[0] * np.cos(obj.vertical_rotation[0]) - obj.vertical_length[1] * np.cos(obj.vertical_rotation[1])
            horizontal_length = np.sqrt(delta_y * delta_y + delta_x * delta_x) + obj.vertical_thickness[0]
            horizontal_rotation = np.arctan(delta_y / delta_x)
            vertical_x_offset = (-obj.vertical_length[0] * np.sin(obj.vertical_rotation[0]) - obj.vertical_length[1] * np.sin(obj.vertical_rotation[1])) / 2
            vertical_y_offset = (obj.vertical_length[1] * np.cos(obj.vertical_rotation[1]) + obj.vertical_length[0] * np.cos(obj.vertical_rotation[0])) / 2
            vertical_mesh_position = [
                vertical_x_offset, 
                vertical_y_offset + obj.horizontal_thickness[0] / 2, 
                0
            ]
            vertical_mesh_rotation = [0, 0, horizontal_rotation]

            affordable_handle = np.array([[-horizontal_length/2, -obj.horizontal_thickness[0]/2, -obj.horizontal_thickness[1]/2],
                                          [horizontal_length/2, obj.horizontal_thickness[0]/2, obj.horizontal_thickness[1]/2]])
            affordable_handle[0] -= AFFORDACE_PROXIMITY_THRES
            affordable_handle[1] += AFFORDACE_PROXIMITY_THRES
            
            _pt = inverse_transformation(pt, obj.position, obj.rotation)
            __pt = inverse_transformation(_pt, vertical_mesh_position, vertical_mesh_rotation)

            return (__pt[0] >= affordable_handle[0, 0] and 
                    __pt[0] <= affordable_handle[1, 0] and 
                    __pt[1] >= affordable_handle[0, 1] and 
                    __pt[1] <= affordable_handle[1, 1] and 
                    __pt[2] >= affordable_handle[0, 2] and 
                    __pt[2] <= affordable_handle[1, 2])
    


        elif (isinstance(obj, Curved_Handle)):
            _pt = inverse_transformation(pt, obj.position, obj.rotation)

            mesh_position = [0, 0, 0]
            mesh_rotation = [-np.pi/2, 0, 0]

            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            _pt_radius = np.linalg.norm(__pt)
            _pt_angle = np.arctan2(__pt[2], __pt[0])

            return (_pt_radius >= obj.radius[1] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.radius[0] + AFFORDACE_PROXIMITY_THRES and
                    _pt_angle >= obj.exist_angle[0]/3 and _pt_angle <= obj.exist_angle[0]*2/3)
        

        elif (isinstance(obj, Round_U_Handle)):
            _pt = inverse_transformation(pt, obj.position, obj.rotation)

            mesh_position = [0, obj.vertical_length[0], 0]
            mesh_rotation = [-np.pi/2, 0, 0]

            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            _pt_radius = np.linalg.norm(__pt)
            _pt_angle = np.arccos(__pt[0] / np.sqrt(__pt[0]**2 + __pt[2]**2))

            return (_pt_radius >= obj.inner_radius[0] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.vertical_separation[0] / 2 + AFFORDACE_PROXIMITY_THRES and
                    _pt_angle >= np.pi/3 and _pt_angle <= np.pi*2/3)
        

        elif (isinstance(obj, Flat_U_Handle)):
            _pt = inverse_transformation(pt, obj.position, obj.rotation)

            mesh_position = [0, obj.vertical_size[1], 0]
            mesh_rotation = [-np.pi/2, 0, 0]

            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            outer_radius = (obj.vertical_separation[0] + obj.vertical_size[0]) / 2
            inner_radius = (obj.vertical_separation[0] - obj.vertical_size[0]) / 2
            _pt_radius = np.linalg.norm(__pt)
            _pt_angle = np.arctan2(__pt[2], __pt[0])

            return (_pt_radius >= inner_radius - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= outer_radius + AFFORDACE_PROXIMITY_THRES and
                    __pt[1] >= -obj.vertical_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.vertical_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_angle >= np.pi/3 and _pt_angle <= np.pi*2/3)
        
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
