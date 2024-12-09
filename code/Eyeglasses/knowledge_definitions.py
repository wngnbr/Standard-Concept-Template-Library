from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def leg_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Regular_Leg)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ')[0]

            leg_interval = obj.offset_x[0] * 2 + obj.glass_interval
            for direction in [-1, 1]:

                position = np.array([direction * leg_interval / 2, 0, 0])
                rotation = np.array([-obj.rotation_1[0], -direction * obj.rotation_1[1], 0])
                __pt = apply_transformation(np.array([_pt]), -position, -rotation, rotation_order='YXZ', offset_first=True)[0]

                position = np.array([0, 0, -obj.size1[2]])
                rotation = np.array([-obj.rotation_2[0], -direction * obj.rotation_2[1], 0])
                __pt = apply_transformation(np.array([__pt]), -position, -rotation, rotation_order='YXZ', offset_first=True)[0]

                mesh_position = np.array([
                    -direction * obj.size2[0] / 2, 
                    -obj.size2[1] / 2, 
                    -obj.size2[2] / 2
                ])
                mesh_rotation = np.array([0, 0, 0])
                __pt = inverse_transformation(__pt, mesh_position, mesh_rotation)

                if (__pt[0] >= -obj.size2[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= obj.size2[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.size2[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.size2[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    __pt[2] >= -obj.size2[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= obj.size2[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
                
            return False
        


        elif (isinstance(obj, Trifold_Leg)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ')[0]

            leg_interval = obj.offset_x[0] * 2 + obj.glass_interval
            for direction in [-1, 1]:

                position = np.array([
                    direction * leg_interval / 2 + obj.offset_x[0] + direction * obj.connector_size[0], 
                    obj.connector_size[1] / 2, 
                    0
                ])
                rotation = np.array([0, 0, 0])
                __pt = apply_transformation(np.array([_pt]), -position, -rotation, offset_first=True)[0]


                position = np.array([0, 0, 0])
                rotation = np.array([-obj.rotation_1[0], -direction * obj.rotation_1[1], 0])
                __pt = apply_transformation(np.array([__pt]), -position, -rotation, offset_first=True)[0]

                position = np.array([0, 0, -obj.size1[2]])
                rotation = np.array([-obj.rotation_2[0], -direction * obj.rotation_2[1], 0])
                __pt = apply_transformation(np.array([__pt]), -position, -rotation, offset_first=True)[0]

                mesh_position = np.array([
                    -direction * obj.size2[0] / 2, 
                    -obj.size2[1] / 2, 
                    -obj.size2[2] / 2
                ])
                mesh_rotation = np.array([0, 0, 0])
                __pt = inverse_transformation(__pt, mesh_position, mesh_rotation)

                if (__pt[0] >= -obj.size2[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= obj.size2[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.size2[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.size2[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    __pt[2] >= -obj.size2[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= obj.size2[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
                
            return False
        
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
