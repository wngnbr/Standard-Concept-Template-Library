from concept_template import *
from geometry_template import *
from knowledge_utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Trifold_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ', offset_first=True)[0]
            
            mesh_position = np.array([0, 0, obj.bottom_size[2] + obj.top_size[2] / 2])
            mesh_rotation = np.array([0, 0, 0])
            
            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            return (__pt[0] >= -obj.top_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[0] <= obj.top_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.top_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.top_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] >= -obj.top_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] <= obj.top_size[2]/2 + AFFORDACE_PROXIMITY_THRES )
    


        elif (isinstance(obj, Claw_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ', offset_first=True)[0]

            for i in range(obj.num_forks[0]):
                rotation_tmp = np.pi * 2 / obj.num_forks[0] * i
                rotate_length = obj.bottom_size[0] + obj.fork_size[0] / 2 * np.cos(obj.fork_tilt_rotation[0])
                mesh_position = np.array([
                    rotate_length * np.cos(rotation_tmp), 
                    rotate_length * np.sin(rotation_tmp), 
                    -obj.fork_size[2] / 2 + obj.bottom_size[1] + obj.fork_size[0] / 2 * np.sin(obj.fork_tilt_rotation[0]) - obj.fork_offset[0]
                ])
                mesh_rotation = np.array([
                    0, 
                    -obj.fork_tilt_rotation[0], 
                    rotation_tmp
                ])
                __pt = apply_transformation(np.array([_pt]), -mesh_position, -mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
                if (__pt[0] >= -obj.fork_size[0]/3 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[0] <= obj.fork_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.fork_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.fork_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] >= -obj.fork_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] <= obj.fork_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
                
            return False
                

        elif (isinstance(obj, Round_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ', offset_first=True)[0]

            circle_mesh_position = np.array([0, 0, obj.bottom_size[1] - obj.fork_size[2] / 2 + obj.fork_size[0] * np.sin(obj.fork_tilt_rotation[0]) - obj.fork_offset[0]])
            circle_mesh_rotation = np.array([np.pi / 2, 0, 0])
            _pt = inverse_transformation(_pt, circle_mesh_position, circle_mesh_rotation)

            outer_radius = obj.circle_size[0]
            inner_radius = obj.bottom_size[0] + obj.fork_size[0] * np.cos(obj.fork_tilt_rotation[0])
            _pt_radius = np.linalg.norm(_pt[[0, 2]])
            
            return (_pt_radius >= inner_radius - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= outer_radius + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.circle_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.circle_size[1]/2 + AFFORDACE_PROXIMITY_THRES)

        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
