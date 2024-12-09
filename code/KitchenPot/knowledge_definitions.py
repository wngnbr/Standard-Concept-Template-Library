from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, Cuboidal_Tophandle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY', offset_first=True)[0]

            mesh_position = np.array([0, obj.size[1] / 2, 0])
            mesh_rotation = np.array([0, 0, 0])
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            return (_pt[0] >= -obj.size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[0] <= obj.size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[1] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    _pt[2] >= -obj.size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[2] <= obj.size[1]/2 + AFFORDACE_PROXIMITY_THRES)



        elif (isinstance(obj, Trifold_Tophandle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZXY', offset_first=True)[0]

            vertical_mesh_position = np.array([
                0, 
                obj.mounting_size[1] * np.cos(obj.mounting_rotation[0]) + obj.grip_size[1] / 2, 
                0
            ])
            vertical_mesh_rotation = np.array([0, 0, 0])
            _pt = inverse_transformation(_pt, vertical_mesh_position, vertical_mesh_rotation)

            return (_pt[0] >= -obj.grip_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[0] <= obj.grip_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.grip_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[1] <= obj.grip_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    _pt[2] >= -obj.grip_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    _pt[2] <= obj.grip_size[1]/2 + AFFORDACE_PROXIMITY_THRES)
        


        elif (isinstance(obj, Semi_Ring_Tophandle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            mesh_position = np.array([
                0, 
                -obj.curve_size[0] * np.cos(obj.curve_exist_angle[0] / 2), 
                0
            ])
            mesh_rotation = np.array([
                -np.pi / 2,
                -np.pi / 2 + obj.curve_exist_angle[0] / 2,
                0
            ])
            _pt = apply_transformation(np.array([_pt]), -mesh_position, -mesh_rotation, rotation_order='ZXY', offset_first=True)[0]

            _pt_radius = np.linalg.norm(_pt[[0, 2]])
            _pt_angle = np.arctan2(_pt[2], _pt[0])

            return (_pt_radius <= obj.curve_size[0] + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius >= obj.curve_size[1] - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= obj.curve_size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] >= -obj.curve_size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt_angle >= obj.curve_exist_angle[0]/3 and _pt_angle <= obj.curve_exist_angle[0]*2/3)
        


        elif (isinstance(obj, Multilevel_Tophandle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), rotation_order='ZYX', offset_first=True)[0]

            delta_height = sum([getattr(obj, 'level_'+ str(i+1) +'_size')[2] for i in range(obj.num_levels[0]-1)]) + \
                getattr(obj, 'level_'+ str(obj.num_levels[0]) +'_size')[2] / 2
            mesh_position = np.array([0, delta_height, 0])
            mesh_rotation = np.array([0, 0, 0])
            _pt = inverse_transformation(_pt, mesh_position, mesh_rotation)

            _pt_radius = np.linalg.norm(_pt[[0, 2]])

            size = getattr(obj, 'level_'+ str(obj.num_levels[0]) +'_size')
            return (_pt[1] >= -size[2]/2 - AFFORDACE_PROXIMITY_THRES and 
                    _pt[1] <= size[2]/2 + AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= (size[0] + size[1]) / 2 - _pt[1] * (size[1] - size[0]) / size[2] + AFFORDACE_PROXIMITY_THRES)

        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
