from concept_template import *
from geometry_template import *
from knowledge_utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Regular_Body)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            base_mesh_position = np.array([0, (obj.top_size[0] - obj.bottom_size[1]) / 2, 0])
            base_mesh_rotation = np.array([0, 0, -np.pi / 2])
            _pt = inverse_transformation(_pt, base_mesh_position, base_mesh_rotation)

            return (
                _pt[0] >= -(obj.top_size[0] + obj.middle_size[0] + obj.bottom_size[1])/3 and
                _pt[0] <= (obj.top_size[0] + obj.middle_size[0] + obj.bottom_size[1])/3 and
                _pt[1] >= -obj.bottom_size[0]/2 and 
                _pt[1] <= obj.bottom_size[0]/2 and 
                (
                    (
                        _pt[2] >= -obj.bottom_size[2]/2 and 
                        _pt[2] <= -(obj.bottom_size[2] - obj.middle_size[1] * 2)/2
                    ) or 
                    (
                        _pt[2] >= (obj.bottom_size[2] - obj.middle_size[1] * 2)/2 and 
                        _pt[2] <= obj.bottom_size[2]/2 
                    )
                )
            )
    
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
