from concept_template import *
from geometry_template import *
from knowledge_utils import *


def cover_affordance(obj, pt):

    def is_affordance(obj, pt):
        if (isinstance(obj, Regular_Cover)):
            _pt = inverse_transformation(pt, obj.position, obj.rotation)

            return (_pt[0] >= -obj.outer_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[0] <= obj.outer_size[0]/2 + AFFORDACE_PROXIMITY_THRES and
                    _pt[2] <= obj.outer_size[2] + AFFORDACE_PROXIMITY_THRES and
                    _pt[2] >= obj.outer_size[2]/2 + obj.inner_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[1] >= -AFFORDACE_PROXIMITY_THRES and
                    _pt[1] <= obj.outer_size[1] + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Fourfold_Cover)):
            _pt = inverse_transformation(pt, obj.position, obj.rotation)


            if (obj.has_cover[0] == 1):                 # front
                cover_position = np.array([
                    0, 
                    obj.front_behind_size[1] * np.cos(obj.cover_rotation[0]) / 2, 
                    obj.cover_separation[0] / 2 + obj.front_behind_size[1] * np.sin(obj.cover_rotation[0]) / 2
                ])
                cover_rotation = np.array([
                    obj.cover_rotation[0],
                    0,
                    0
                ])
                __pt = inverse_transformation(_pt, cover_position, cover_rotation)

                if (__pt[0] >= -obj.front_behind_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= obj.front_behind_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= obj.front_behind_size[1]/8 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.front_behind_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    __pt[2] >= -obj.front_behind_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= obj.front_behind_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
            

            if (obj.has_cover[1] == 1):               # back
                cover_position = np.array([
                    0, 
                    obj.front_behind_size[1] * np.cos(obj.cover_rotation[1]) / 2, 
                    -obj.cover_separation[0] / 2 + obj.front_behind_size[1] * np.sin(obj.cover_rotation[1]) / 2
                ])
                cover_rotation = np.array([
                    obj.cover_rotation[1],
                    0,
                    0
                ])
                __pt = inverse_transformation(_pt, cover_position, cover_rotation)

                if (__pt[0] >= -obj.front_behind_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= obj.front_behind_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= obj.front_behind_size[1]/8 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.front_behind_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    __pt[2] >= -obj.front_behind_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= obj.front_behind_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
                

            if (obj.has_cover[2] == 1):               # left
                cover_position = np.array([
                    -obj.cover_separation[1] / 2 - obj.left_right_size[1] * np.sin(obj.cover_rotation[2]) / 2, 
                    obj.left_right_size[1] * np.cos(obj.cover_rotation[2]) / 2, 
                    0
                ])
                cover_rotation = np.array([
                    0, 
                    0,
                    obj.cover_rotation[2],
                ])
                __pt = inverse_transformation(_pt, cover_position, cover_rotation)

                if (__pt[0] >= -obj.left_right_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= obj.left_right_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= obj.left_right_size[1]/8 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.left_right_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    __pt[2] >= -obj.left_right_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= obj.left_right_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
                

            if (obj.has_cover[3] == 1):               # right
                cover_position = np.array([
                    obj.cover_separation[1] / 2 - obj.left_right_size[1] * np.sin(obj.cover_rotation[3]) / 2, 
                    obj.left_right_size[1] * np.cos(obj.cover_rotation[3]) / 2, 
                    0
                ])
                cover_rotation = np.array([
                    0, 
                    0,
                    obj.cover_rotation[3],
                ])
                __pt = inverse_transformation(_pt, cover_position, cover_rotation)

                if (__pt[0] >= -obj.left_right_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= obj.left_right_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= obj.left_right_size[1]/8 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.left_right_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    __pt[2] >= -obj.left_right_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= obj.left_right_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
                
            return False

        return False
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
