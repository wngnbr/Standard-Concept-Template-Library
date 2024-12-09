from concept_template import *
from geometry_template import *
from knowledge_utils import *


def cover_affordance(obj, pt):

    def is_affordance(obj, pt):
        if (isinstance(obj, Cylindrical_Cover)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            return (_pt_radius <= obj.size[0] + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Cuboidal_Cover)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = [0, obj.height[0] / 2, obj.bottom_size[1] / 2]
            mesh_rotation = [0, 0, 0]
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            return (_pt[1] >= -obj.height[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.height[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[0] <= ((obj.bottom_size[0] + obj.top_size[0]) / 2 - _pt[1] * (obj.bottom_size[0] - obj.top_size[0]) / obj.height[0]) / 2 + obj.top_offset[0] * (_pt[1] + obj.height[0]/2) / obj.height[0] + AFFORDACE_PROXIMITY_THRES and 
                    _pt[0] >= -((obj.bottom_size[0] + obj.top_size[0]) / 2 - _pt[1] * (obj.bottom_size[0] - obj.top_size[0]) / obj.height[0]) / 2 + obj.top_offset[0] * (_pt[1] + obj.height[0]/2) / obj.height[0] - AFFORDACE_PROXIMITY_THRES and 
                    _pt[2] <= ((obj.bottom_size[1] + obj.top_size[1]) / 2 - _pt[1] * (obj.bottom_size[1] - obj.top_size[1]) / obj.height[0]) / 2 + obj.top_offset[1] * (_pt[1] + obj.height[0]/2) / obj.height[0] + AFFORDACE_PROXIMITY_THRES and 
                    _pt[2] >= ((obj.bottom_size[1] + obj.top_size[1]) / 2 - _pt[1] * (obj.bottom_size[1] - obj.top_size[1]) / obj.height[0]) * 2 / 5 + obj.top_offset[1] * (_pt[1] + obj.height[0]/2) / obj.height[0] - AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Double_Layer_Cuboidal_Cover)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            bottom_mesh_position = [
                0,
                obj.bottom_size[1] / 2,
                obj.bottom_size[2] / 2
            ]
            bottom_mesh_rotation = [0, 0, 0]
            _pt = inverse_transformation(_pt, bottom_mesh_position, bottom_mesh_rotation)

            return (_pt[0] >= -obj.bottom_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[0] <= obj.bottom_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.bottom_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.bottom_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[2] >= max(obj.bottom_size[2]*2/5, obj.top_size[2]/2) - AFFORDACE_PROXIMITY_THRES and 
                    _pt[2] <= obj.bottom_size[2]/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Cylindrical_Hollow_Cover)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = [0, obj.size[2] / 2, 0]
            mesh_rotation = [0, 0, 0]
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            return (_pt[1] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius >= obj.size[1] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.size[0] + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Cuboidal_Hollow_Cover)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = [0, obj.outer_size[1] / 2, 0]
            mesh_rotation = [0, 0, 0]
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            return (_pt[1] >= -obj.outer_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.outer_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                    ((
                        (_pt[0] >= -obj.outer_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                         _pt[0] <= -obj.inner_size[0]/2 + obj.inner_offset[0] + AFFORDACE_PROXIMITY_THRES) 
                        or
                        (_pt[0] >= obj.inner_size[0]/2 + obj.inner_offset[0] - AFFORDACE_PROXIMITY_THRES and
                         _pt[0] <= obj.outer_size[0]/2 + AFFORDACE_PROXIMITY_THRES) 
                    ) or 
                    (
                        (_pt[2] >= -obj.outer_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                         _pt[2] <= -obj.inner_size[1]/2 + obj.inner_offset[1] + AFFORDACE_PROXIMITY_THRES) 
                        or
                        (_pt[2] >= obj.inner_size[1]/2 + obj.inner_offset[1] - AFFORDACE_PROXIMITY_THRES and
                         _pt[2] <= obj.outer_size[2]/2 + AFFORDACE_PROXIMITY_THRES) 
                    ) 
                    ))
        

        return False
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
