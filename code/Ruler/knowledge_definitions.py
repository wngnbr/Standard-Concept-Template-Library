from concept_template import *
from geometry_template import *
from knowledge_utils import *


def ruler_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Symmetrical_body)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            mesh_1_position = np.array([
                -obj.separation[0] - (obj.size[0] * np.cos(obj.body_rotation[0]) + obj.size[1] * np.sin(obj.body_rotation[0])) / 2,
                (obj.size[1] * np.cos(obj.body_rotation[0]) - obj.size[0] * np.sin(obj.body_rotation[0])) / 2, 
                -obj.left_right_offset[0]
            ])
            mesh_1_rotation = np.array([0, 0, obj.body_rotation[0]])

            __pt = inverse_transformation(_pt, mesh_1_position, mesh_1_rotation)

            if (__pt[0] >= -obj.size[0]*2/5 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= obj.size[0]*2/5 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >=  - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            

            mesh_2_position = np.array([
                obj.separation[0] + (obj.size[0] * np.cos(obj.body_rotation[0]) + obj.size[1] * np.sin(obj.body_rotation[0])) / 2,
                (obj.size[1] * np.cos(obj.body_rotation[0]) - obj.size[0] * np.sin(obj.body_rotation[0])) / 2, 
                -obj.left_right_offset[0]
            ])
            mesh_2_rotation = np.array([0, 0, -obj.body_rotation[0]])

            __pt = inverse_transformation(_pt, mesh_2_position, mesh_2_rotation)

            if (__pt[0] >= -obj.size[0]*2/5 - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= obj.size[0]*2/5 + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >=  - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            
            return False
        

        elif (isinstance(obj, Asymmetrical_body)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]
            
            shorter_size = obj.left_size if obj.left_size[0] < obj.right_size[0] else obj.right_size

            mesh_1_position = np.array([
                -obj.separation[0] - obj.left_size[0] / 2 * np.cos(obj.body_rotation[0]) - obj.left_size[1] / 2 * np.sin(obj.body_rotation[0]),
                -obj.left_size[0] / 2 * np.sin(obj.body_rotation[0]) + obj.left_size[1] / 2 * np.cos(obj.body_rotation[0]), 
                -obj.left_right_offset[0]
            ])
            mesh_1_rotation = np.array([0, 0, obj.body_rotation[0]])

            __pt = inverse_transformation(_pt, mesh_1_position, mesh_1_rotation)

            if (__pt[0] >= -(obj.left_size[0]/2 - shorter_size[0]/10) - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= (obj.left_size[0]/2 - shorter_size[0]/10) + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >=  - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.left_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.left_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.left_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            

            mesh_2_position = np.array([
                obj.separation[0] + obj.right_size[0] / 2 * np.cos(obj.body_rotation[0]) - obj.right_size[1] / 2 * np.sin(-obj.body_rotation[0]),
                obj.right_size[0] / 2 * np.sin(-obj.body_rotation[0]) + obj.right_size[1] / 2 * np.cos(obj.body_rotation[0]), 
                obj.left_right_offset[0]
            ])
            mesh_2_rotation = np.array([0, 0, -obj.body_rotation[0]])

            __pt = inverse_transformation(_pt, mesh_2_position, mesh_2_rotation)

            if (__pt[0] >= -(obj.right_size[0]/2 - shorter_size[0]/10) - AFFORDACE_PROXIMITY_THRES and 
                __pt[0] <= (obj.right_size[0]/2 - shorter_size[0]/10) + AFFORDACE_PROXIMITY_THRES and 
                __pt[1] >=  - AFFORDACE_PROXIMITY_THRES and 
                __pt[1] <= obj.right_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                __pt[2] >= -obj.right_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                __pt[2] <= obj.right_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                return True
            
            return False
    
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
