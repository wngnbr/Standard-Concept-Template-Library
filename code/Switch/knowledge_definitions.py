from concept_template import *
from geometry_template import *
from knowledge_utils import *


def switch_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Round_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            for i in range(obj.number_of_switch[0]):
                mesh_position = [
                    getattr(obj, 'offset_%d'%(i+1))[0],
                    getattr(obj, 'offset_%d'%(i+1))[1],
                    obj.offset_Z[0]
                ]
                mesh_rotation = [np.pi / 2 + obj.switch_rotation[0], 0, 0]
                __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)
                _pt_radius = np.linalg.norm(__pt[[0, 2]])

                if (_pt_radius <= obj.size[0] + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
            
            return False
        

        elif (isinstance(obj, FlipX_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            for i in range(obj.number_of_switch[0]):
                mesh_position = [(obj.separation[0] + obj.size[0]) * i, 0, 0]
                mesh_rotation = [obj.switch_rotation[0], 0, 0]
                __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

                if (obj.switch_rotation[0] >= 0):
                    if (__pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                        __pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= - AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                        __pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                        return True
                
                elif (obj.switch_rotation[0] < 0):
                    if (__pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                        __pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] <= AFFORDACE_PROXIMITY_THRES and 
                        __pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                        __pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                        return True
                
            return False
        

        elif (isinstance(obj, FlipY_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            for i in range(obj.number_of_switch[0]):
                mesh_position = [(obj.separation[0] + obj.size[0]) * i, 0, 0]
                mesh_rotation = [0, obj.switch_rotation[0], 0]
                __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

                if (obj.switch_rotation[0] >= 0):
                    if (__pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and 
                        __pt[0] <= AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                        __pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                        return True
                
                elif (obj.switch_rotation[0] < 0):
                    if (__pt[0] >= - AFFORDACE_PROXIMITY_THRES and 
                        __pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                        __pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                        return True
                
            return False
        

        elif (isinstance(obj, Lever_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            for i in range(obj.number_of_switch[0]):
                mesh_position = [(obj.separation[0] + obj.main_size[0]) * i, 0, obj.main_size[2] / 2,]
                mesh_rotation = [np.pi / 2, 0, 0]

                position = np.array([
                    obj.inter_offset[0],
                    obj.inter_offset[1],
                    obj.inter_offset[2] + obj.base_size[1],
                ])
                rotation = np.array([0, 0, 0])
                __pt = apply_transformation(np.array([_pt]), -position, -rotation, rotation_order='ZYX', offset_first='True')[0]
                position = np.array([0, 0, 0])
                rotation = np.array([obj.switch_rotation[0], 0, 0])
                __pt = apply_transformation(np.array([__pt]), -position, -rotation, rotation_order='ZYX', offset_first='True')[0]
                __pt = inverse_transformation(__pt, mesh_position, mesh_rotation)

                _pt_radius = np.linalg.norm(__pt[[0, 2]])

                if (__pt[1] >= - AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] <= obj.main_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= (obj.main_size[0] + obj.main_size[1]) / 2 - _pt[1] * (obj.main_size[1] - obj.main_size[0]) / obj.main_size[2] + AFFORDACE_PROXIMITY_THRES):
                    return True
                
            return False
        
        
        else:
            return False
    
    
    return is_affordance(obj, pt)



def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
