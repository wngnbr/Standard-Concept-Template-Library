from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def sphere_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Standard_Sphere)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            horizontal_angle = np.arccos(_pt[0] / obj.radius[0])
            tilt_angle = np.arctan2(_pt[1], obj.radius[0])

            return (horizontal_angle >= 0 and
                    horizontal_angle <= np.pi/6 and 
                    tilt_angle >= -np.pi/6 and 
                    tilt_angle <= np.pi/6)
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
