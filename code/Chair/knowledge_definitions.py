from concept_template import *
from geometry_template import *
from knowledge_utils import *


def seat_affordance(obj, pt):

    def is_affordance(obj, pt):
        
        if (isinstance(obj, Regular_seat)):
            _pt = inverse_transformation(pt, obj.position, obj.rotation)
            
            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            return (_pt_radius <= min(obj.size[0], obj.size[2])/3 and
                    _pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES)
    


        elif (isinstance(obj, Round_seat)):
            _pt = inverse_transformation(pt, obj.position, obj.rotation)
            
            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            return (_pt_radius <= obj.radius*2/3 and
                    _pt[1] >= -obj.height/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[1] <= obj.height/2 + AFFORDACE_PROXIMITY_THRES)
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
