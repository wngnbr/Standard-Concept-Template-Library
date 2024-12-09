from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def nozzle_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Press_Nozzle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ', offset_first=True)[0]

            delta_height = sum([getattr(obj, 'level_' + str(i+1) + '_size')[1] for i in range(obj.num_levels[0]-1)])
            delta_height += getattr(obj, 'level_' + str(obj.num_levels[0]) + '_size')[1] / 2
            mesh_position = np.array([0, delta_height, 0])
            mesh_rotation = np.array([0, 0, 0])

            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            _pt_radius = np.linalg.norm(__pt[[0, 2]])

            return (_pt_radius <= getattr(obj, 'level_' + str(obj.num_levels[0]) + '_size')[0] and
                    __pt[1] >= getattr(obj, 'level_' + str(obj.num_levels[0]) + '_size')[1]*3/10 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= getattr(obj, 'level_' + str(obj.num_levels[0]) + '_size')[1]/2 + AFFORDACE_PROXIMITY_THRES)
        

            # total_offset_z = obj.nozzle_length[0] / 2 * np.cos(obj.nozzle_rotation[0]) + getattr(obj, 'level_%d_size'%(obj.num_levels[0]))[0]
            # mesh_position = np.array([
            #     0, 
            #     delta_height + obj.nozzle_offset[0] - obj.nozzle_length[0] / 2 * np.sin(obj.nozzle_rotation[0]),  
            #     total_offset_z
            # ])
            # mesh_rotation = np.array([obj.nozzle_rotation[0], 0, 0])

            # __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            # return (__pt[0] >= -obj.nozzle_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
            #         __pt[0] <= obj.nozzle_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
            #         __pt[1] >= -obj.nozzle_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
            #         __pt[1] <= obj.nozzle_size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
            #         __pt[2] >= -obj.nozzle_length[0]/2 - AFFORDACE_PROXIMITY_THRES and 
            #         __pt[2] <= obj.nozzle_length[0]/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Spray_Nozzle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ', offset_first=True)[0]
            
            handle_total_offset_z = obj.top_offset[0] + obj.top_size[2] / 2 + obj.handle_offset[0] - \
                obj.handle_size[2] / 2 - obj.handle_size[1] / 2 * np.sin(obj.handle_rotation[0])
            handle_total_offset_y = -obj.handle_size[1] / 2 * np.cos(obj.handle_rotation[0]) - obj.top_size[1] / 2
            handle_mesh_position_1 = np.array([
                0, 
                handle_total_offset_y,
                handle_total_offset_z
            ])
            handle_mesh_rotation = np.array([obj.top_rotation[0], 0, 0])
            handle_mesh_position_1 = adjust_position_from_rotation(handle_mesh_position_1, handle_mesh_rotation)
            delta_height = obj.bottom_size[1] + obj.middle_size[1] + obj.top_size[1]/2
            top_offset_y = -obj.top_offset[0] * np.tan(obj.top_rotation[0])
            handle_mesh_position_2 = np.array([
                0, 
                delta_height + top_offset_y,
                0
            ])
            handle_mesh_position = handle_mesh_position_1 + handle_mesh_position_2
            handle_mesh_rotation[0] += obj.handle_rotation[0]
            
            __pt = inverse_transformation(_pt, handle_mesh_position, handle_mesh_rotation)

            return (__pt[0] >= -obj.handle_size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[0] <= obj.handle_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.handle_size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= 0 and 
                    __pt[2] >= -obj.handle_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[2] <= obj.handle_size[2]/2 + AFFORDACE_PROXIMITY_THRES )
    
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
