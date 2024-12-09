from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def support_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, CuboidalRear_Support)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ')[0]

            for i in range(obj.number_of_supports[0]):
                mesh_position = np.array([
                    (obj.separation[0] + obj.size[0]) * i, 
                    -obj.size[1] / 2, 
                    -obj.size[2] / 2
                ])
                mesh_rotation = np.array([0, 0, 0])

                __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

                if (__pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= -obj.size[1]/4 + AFFORDACE_PROXIMITY_THRES and
                    __pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True

            return False
        


        elif (isinstance(obj, Cuboidal_Support)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ')[0]

            for i in range(obj.number_of_supports[0]):
                mesh_position = np.array([
                    (obj.separation[0] + obj.size[0]) * i, 
                    -obj.size[1] / 2,
                    0
                ])
                mesh_rotation = np.array([0, 0, 0])

                __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

                if (__pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= AFFORDACE_PROXIMITY_THRES and
                    __pt[2] >= -obj.size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= obj.size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True

            return False
        

        elif (isinstance(obj, Trifold_Support)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ')[0]


            if obj.has_middle_part[0]:
                base_position = [0, 0, 0]
                base_position = [
                    base_position[0] + obj.upper_offset[0],
                    base_position[1] + obj.upper_offset[1],
                    base_position[2] + obj.upper_offset[2] - obj.upper_part_size[2],
                ]
                mesh_position = np.array([
                    base_position[0],
                    base_position[1] + obj.middle_offset[0] - obj.middle_part_size[1] / 2,
                    base_position[2] - obj.middle_part_size[2] / 2
                ])
                mesh_rotation = np.array([0, 0, 0])

                __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

                return (__pt[0] >= -obj.middle_part_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[0] <= obj.middle_part_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= -obj.middle_part_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[1] <= obj.middle_part_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                        __pt[2] >= -obj.middle_part_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[2] <= obj.middle_part_size[2]/2 + AFFORDACE_PROXIMITY_THRES)
            

            if obj.has_bottom_part[0]:
                base_position = [0, 0, 0]
                base_position = [
                    base_position[0] + obj.upper_offset[0],
                    base_position[1] + obj.upper_offset[1],
                    base_position[2] + obj.upper_offset[2] - obj.upper_part_size[2],
                ]
                base_position = [
                    base_position[0],
                    base_position[1] + obj.middle_offset[0] - obj.middle_part_size[1],
                    base_position[2] - obj.middle_part_size[2],
                ]
                base_position = np.array([
                    base_position[0],
                    base_position[1],
                    base_position[2] + obj.bottom_offset[0],
                ])
                mesh_position = np.array([
                    0, 
                    -obj.bottom_part_size[1] / 2, 
                    obj.bottom_part_size[2] / 2
                ])
                
                mesh_rotation_1 = np.array([obj.rotation_of_bottom[0], 0, 0])
                mesh_rotation_2 = np.array([0, 0, 0])

                __pt = apply_transformation(_pt, -base_position, -mesh_rotation_1, offset_first=True)
                __pt = inverse_transformation(__pt, mesh_position, mesh_rotation_2)

                return (__pt[0] >= -obj.bottom_part_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[0] <= obj.bottom_part_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= -obj.bottom_part_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[1] <= obj.bottom_part_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                        __pt[2] >= -obj.bottom_part_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[2] <= obj.bottom_part_size[2]/2 + AFFORDACE_PROXIMITY_THRES)


            if obj.has_upper_part[0]:
                mesh_position = np.array([
                    obj.upper_offset[0],
                    obj.upper_offset[1],
                    obj.upper_offset[2] - obj.upper_part_size[2] / 2
                ])
                mesh_rotation = np.array([0, 0, 0])

                __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

                return (__pt[0] >= -obj.upper_part_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[0] <= obj.upper_part_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= -obj.upper_part_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[1] <= obj.upper_part_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                        __pt[2] >= -obj.upper_part_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[2] <= obj.upper_part_size[2]/2 + AFFORDACE_PROXIMITY_THRES)


            return False
        

        elif (isinstance(obj, TShaped_Support)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='YXZ')[0]

            mesh_position = np.array([
                obj.upper_offset[0],
                obj.upper_offset[1] - obj.middle_part_size[1] / 2 - obj.upper_part_size[1] / 2,
                obj.upper_offset[2] - obj.middle_part_size[2] / 2 + obj.middle_offset[0],
            ])
            mesh_rotation = np.array([0, 0, 0])
            __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            return (__pt[0] >= -obj.middle_part_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= obj.middle_part_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.middle_part_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.middle_part_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    __pt[2] >= -obj.middle_part_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= obj.middle_part_size[2]/2 + AFFORDACE_PROXIMITY_THRES)
    
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
