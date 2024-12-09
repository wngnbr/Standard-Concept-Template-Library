from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def switch_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, SimplifiedZ_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            for offset in obj.offsets:
                
                mesh_position = np.array(offset)
                mesh_rotation = np.array([np.pi / 2, 0, 0])
                __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

                _pt_radius = np.linalg.norm(__pt[[0, 2]])

                if (_pt_radius <= obj.size[0] + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True
                
            return False
        


        elif (isinstance(obj, Knob_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            if (len(obj.offsets) > 1):
                for offset in obj.offsets:
                    mesh_position = np.array(offset[-1])
                    mesh_rotation = np.array([0, 0, 0])
                    __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

                    _pt_radius = np.linalg.norm(__pt[[0, 2]])

                    if (_pt_radius <= obj.sizes[-1][0] + AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= - AFFORDACE_PROXIMITY_THRES and
                        __pt[1] <= obj.sizes[-1][1]/2 + AFFORDACE_PROXIMITY_THRES):
                        return True
            
            else:
                for offset in obj.offsets:
                    for i in range(obj.number_of_cylinder[0]):
                        mesh_position = np.array(offset[i])
                        mesh_rotation = np.array([0, 0, 0])
                        __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

                        _pt_radius = np.linalg.norm(__pt[[0, 2]])

                        if (_pt_radius <= obj.sizes[i][0] + AFFORDACE_PROXIMITY_THRES and 
                            __pt[1] >= - AFFORDACE_PROXIMITY_THRES and
                            __pt[1] <= obj.sizes[i][1]/2 + AFFORDACE_PROXIMITY_THRES):
                            return True
                
            return False
        

        elif (isinstance(obj, HandleY_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            for offset in obj.offsets:

                mesh_position = np.array(offset[-1])
                mesh_rotation = np.array([0, 0, 0])
                __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

                if (offset[-1][0] >= 0):
                    if (__pt[0] >= - AFFORDACE_PROXIMITY_THRES and
                        __pt[0] <= obj.sizes[-1][0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= -obj.sizes[-1][1]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[1] <= obj.sizes[-1][1]/2 + AFFORDACE_PROXIMITY_THRES and
                        __pt[2] >= -obj.sizes[-1][2]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[2] <= obj.sizes[-1][2]/2 + AFFORDACE_PROXIMITY_THRES):
                        return True
                else:
                    if (__pt[0] >= -obj.sizes[-1][0]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[0] <= AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= -obj.sizes[-1][1]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[1] <= obj.sizes[-1][1]/2 + AFFORDACE_PROXIMITY_THRES and
                        __pt[2] >= -obj.sizes[-1][2]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[2] <= obj.sizes[-1][2]/2 + AFFORDACE_PROXIMITY_THRES):
                        return True
                
            return False
        

        elif (isinstance(obj, HandleZ_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            for offset in obj.offsets:

                mesh_position = np.array(offset[-1])
                mesh_rotation = np.array([0, 0, 0])
                __pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

                if (offset[-1][0] >= 0):
                    if (__pt[0] >= - AFFORDACE_PROXIMITY_THRES and
                        __pt[0] <= obj.sizes[-1][0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= -obj.sizes[-1][1]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[1] <= obj.sizes[-1][1]/2 + AFFORDACE_PROXIMITY_THRES and
                        __pt[2] >= -obj.sizes[-1][2]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[2] <= obj.sizes[-1][2]/2 + AFFORDACE_PROXIMITY_THRES):
                        return True
                else:
                    if (__pt[0] >= -obj.sizes[-1][0]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[0] <= AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= -obj.sizes[-1][1]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[1] <= obj.sizes[-1][1]/2 + AFFORDACE_PROXIMITY_THRES and
                        __pt[2] >= -obj.sizes[-1][2]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[2] <= obj.sizes[-1][2]/2 + AFFORDACE_PROXIMITY_THRES):
                        return True
                
            return False
        

        elif (isinstance(obj, RegularY_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            sub_mesh_position = np.array([0, obj.size[1] / 2, 0])
            sub_mesh_rotation = np.array([obj.rotation_X[0], obj.rotation_Y[0], 0])
            __pt = apply_transformation(np.array([_pt]), -sub_mesh_position, -sub_mesh_rotation, offset_first=True)[0]

            sub_mesh_position = np.array([0, obj.sub_offset[0], (obj.size[0] ** 2 - obj.sub_size[0] ** 2) ** 0.5])
            sub_mesh_rotation = np.array([obj.sub_rotation[0], 0, 0])
            __pt = apply_transformation(np.array([__pt]), -sub_mesh_position, -sub_mesh_rotation, offset_first=True)[0]

            sub_mesh_position = np.array([0, 0, obj.sub_size[1] / 2 + obj.sub_offset[1]])
            sub_mesh_rotation = np.array([np.pi/2, 0, 0])
            __pt = inverse_transformation(__pt, sub_mesh_position, sub_mesh_rotation)

            _pt_radius = np.linalg.norm(__pt[[0, 2]])

            return (_pt_radius <= obj.sub_size[0] + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.sub_size[1]/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, RegularX_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            sub_mesh_position = np.array([obj.size[1] / 2, 0, 0])
            sub_mesh_rotation = np.array([obj.rotation_X[0], 0, obj.rotation_Z[0]])
            __pt = apply_transformation(np.array([_pt]), -sub_mesh_position, -sub_mesh_rotation, offset_first=True)[0]

            sub_mesh_position = np.array([obj.sub_offset[0], (obj.size[0] ** 2 - obj.sub_size[0] ** 2) ** 0.5, 0])
            sub_mesh_rotation = np.array([0, 0, obj.sub_rotation[0]])
            __pt = apply_transformation(np.array([__pt]), -sub_mesh_position, -sub_mesh_rotation, offset_first=True)[0]

            sub_mesh_position = np.array([0, obj.sub_size[1] / 2 + obj.sub_offset[1], 0])
            sub_mesh_rotation = np.array([0, 0, 0])
            __pt = inverse_transformation(__pt, sub_mesh_position, sub_mesh_rotation)

            _pt_radius = np.linalg.norm(__pt[[0, 2]])

            return (_pt_radius <= obj.sub_size[0] + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.sub_size[1]/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, RegularZ_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            sub_mesh_position = np.array([0, 0, obj.size[1] / 2])
            sub_mesh_rotation = np.array([obj.rotation_X[0], 0, obj.rotation_Z[0]])
            __pt = apply_transformation(np.array([_pt]), -sub_mesh_position, -sub_mesh_rotation, offset_first=True)[0]

            sub_mesh_position = np.array([0, (obj.size[0] ** 2 - obj.sub_size[0] ** 2) ** 0.5, obj.sub_offset[0]])
            sub_mesh_rotation = np.array([obj.sub_rotation[0], 0, 0])
            __pt = apply_transformation(np.array([__pt]), -sub_mesh_position, -sub_mesh_rotation, offset_first=True)[0]

            sub_mesh_position = np.array([0, obj.sub_size[1] / 2 + obj.sub_offset[1], 0])
            sub_mesh_rotation = np.array([0, 0, 0])
            __pt = inverse_transformation(__pt, sub_mesh_position, sub_mesh_rotation)

            _pt_radius = np.linalg.norm(__pt[[0, 2]])

            return (_pt_radius <= obj.sub_size[0] + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.sub_size[1]/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, TShaped_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            mesh_position = np.array([0, 0, 0])
            mesh_rotation = np.array([0, obj.switch_rotation[0], 0])
            _pt = apply_transformation(np.array([_pt]), -mesh_position, -mesh_rotation, offset_first=True)[0]

            mesh_position = np.array(obj.offsets[2])
            mesh_rotation = np.array([0, 0, 0])
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            return (_pt[0] >= -obj.size_3[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[0] <= obj.size_3[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.size_3[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[1] <= obj.size_3[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    _pt[2] >= - AFFORDACE_PROXIMITY_THRES and
                    _pt[2] <= obj.size_3[2]/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Lever_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            start_ends = [
                [obj.positions[0], obj.positions[1]],
                [obj.positions[1], obj.positions[2]],
            ]
            vector = np.array(start_ends[-1][1]) - np.array(start_ends[-1][0])
            for i in range(3):
                if vector[i] <= 0.00001:
                    vector[i] += 0.00001
            length = np.linalg.norm(vector)

            mesh_position = np.array([
                start_ends[-1][0][0],
                start_ends[-1][0][1] + obj.size[1], 
                start_ends[-1][0][2] + obj.size[0]
            ])
            mesh_rotation = np.array([np.arccos(vector[1] / length), np.arctan(vector[0] / vector[2]), 0])

            __pt = apply_transformation(np.array([_pt]), -mesh_position, -mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
            mesh_position = np.array([0, length / 2, 0])
            mesh_rotation = np.array([0, 0, 0])
            __pt = inverse_transformation(__pt, mesh_position, mesh_rotation)

            _pt_radius = np.linalg.norm(__pt[[0, 2]])

            return (_pt_radius <= obj.R[0] + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= length/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, Cuboidal_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            sub_mesh_position = np.array([0, obj.size[1] / 2, 0])
            sub_mesh_rotation = np.array([0, obj.switch_rotation[0], 0])
            _pt = apply_transformation(np.array([_pt]), -sub_mesh_position, -sub_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
            offset_z = obj.sub_size[2] / 2 + obj.sub_offset[1]
            sub_mesh_position = np.array([
                0, 
                offset_z * np.sin(obj.sub_rotation[0]) + obj.sub_offset[0] + obj.size[1] / 2 - obj.sub_size[1] / 2 * np.cos(obj.sub_rotation[0]), 
                offset_z * np.cos(obj.sub_rotation[0])
            ])
            sub_mesh_rotation = np.array([-obj.sub_rotation[0], 0, 0])
            _pt = inverse_transformation(_pt, sub_mesh_position, sub_mesh_rotation)

            return (_pt[0] >= -obj.sub_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[0] <= obj.sub_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.sub_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[1] <= obj.sub_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    _pt[2] >= - AFFORDACE_PROXIMITY_THRES and
                    _pt[2] <= obj.sub_size[2]/2 + AFFORDACE_PROXIMITY_THRES)
        

        elif (isinstance(obj, RotaryX_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            rotations = [obj.rotation0[0], obj.rotation1[0]]

            j = 0
            for existence in [obj.existence_of_switch[0] * -1, obj.existence_of_switch[1]]:
                if existence == 0:
                    j += 1
                    continue

                for i in range(obj.number_of_sub[0]):
                    sub_mesh_position = np.array([
                        obj.offset_x[j] + existence * (obj.main_size_1[1] + obj.main_size_2[1] - obj.sub_size[0] / 2 + obj.sub_offset[0] + obj.interval / 2),
                        0,
                        0,
                    ])
                    sub_mesh_rotation = np.array([rotations[j] + np.pi * 2 * i / obj.number_of_sub[0], 0, 0])
                    __pt = apply_transformation(np.array([_pt]), -sub_mesh_position, -sub_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
                    sub_mesh_position = np.array([0, np.cos(obj.tilt_angle[0]) * obj.sub_size[1] / 2, 0])
                    sub_mesh_rotation = np.array([0, 0, -obj.tilt_angle[0]])
                    __pt = inverse_transformation(__pt, sub_mesh_position, sub_mesh_rotation)

                    if (__pt[0] >= -obj.sub_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[0] <= obj.sub_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= - AFFORDACE_PROXIMITY_THRES and
                        __pt[1] <= obj.sub_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                        __pt[2] >= -obj.sub_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[2] <= obj.sub_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                        return True
                
                j += 1
                    
            return False
        


        elif (isinstance(obj, RotaryY_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            rotations = [obj.rotation0[0], obj.rotation1[0]]

            for j in range(obj.number_of_switch[0]):
                for i in range(obj.number_of_sub[0]):
                    sub_mesh_position = np.array([
                        obj.offset_x[j], 
                        obj.main_size_1[1] + obj.main_size_2[1] - obj.sub_size[1] / 2 + obj.sub_offset[0],
                        0,
                    ])
                    sub_mesh_rotation = np.array([0, rotations[j] + np.pi * 2 * i / obj.number_of_sub[0], 0])
                    __pt = apply_transformation(np.array([_pt]), -sub_mesh_position, -sub_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
                    sub_mesh_position = np.array([0, 0, np.cos(obj.tilt_angle[0]) * obj.sub_size[2] / 2])
                    sub_mesh_rotation = np.array([-obj.tilt_angle[0], 0, 0])
                    __pt = inverse_transformation(__pt, sub_mesh_position, sub_mesh_rotation)

                    if (__pt[0] >= -obj.sub_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[0] <= obj.sub_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= -obj.sub_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[1] <= obj.sub_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                        __pt[2] >= - AFFORDACE_PROXIMITY_THRES and
                        __pt[2] <= obj.sub_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                        return True
                    
            return False
        


        elif (isinstance(obj, RotaryZ_Switch)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY')[0]

            rotations = [obj.rotation0[0], obj.rotation1[0]]

            for j in range(obj.number_of_switch[0]):
                for i in range(obj.number_of_sub[0]):
                    sub_mesh_position = np.array([
                        obj.offset_x[j], 
                        0, 
                        obj.main_size_1[1] + obj.main_size_2[1] - obj.sub_size[2] / 2 + obj.sub_offset[0],
                    ])
                    sub_mesh_rotation = np.array([0, 0, rotations[j] + np.pi * 2 * i / obj.number_of_sub[0]])
                    __pt = apply_transformation(np.array([_pt]), -sub_mesh_position, -sub_mesh_rotation, rotation_order='ZYX', offset_first=True)[0]
            
                    sub_mesh_position = np.array([np.cos(obj.tilt_angle[0]) * obj.sub_size[0] / 2, 0, 0])
                    sub_mesh_rotation = np.array([0, -obj.tilt_angle[0], 0])
                    __pt = inverse_transformation(__pt, sub_mesh_position, sub_mesh_rotation)

                    if (__pt[0] >= - AFFORDACE_PROXIMITY_THRES and
                        __pt[0] <= obj.sub_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                        __pt[1] >= -obj.sub_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[1] <= obj.sub_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                        __pt[2] >= -obj.sub_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                        __pt[2] <= obj.sub_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                        return True
                    
            return False
        
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
