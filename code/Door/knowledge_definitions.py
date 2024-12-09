from concept_template import *
from geometry_template import *
from knowledge_utils import *
from utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):

        if (isinstance(obj, LShape_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), offset_first=True)[0]

            for direction_setting in obj.direction_settings:
                __pt = apply_transformation(
                    np.array([_pt]), 
                    -np.array([direction_setting["handle_y_axis"], 0, 0]),
                    -np.array([0, direction_setting["handle_y_rotation"], 0]),
                    offset_first=True
                )[0]
                __pt = apply_transformation(
                    np.array([__pt]), 
                    -np.array([-direction_setting["handle_y_axis"], 0, 0]),
                    np.array([0, 0, 0]),
                    offset_first=True
                )[0]

                top_mesh_position = np.array([
                    direction_setting["handle_x_position"]
                    + obj.interpiece_offset[0]
                    + direction_setting["handle_x_direction"]
                    * (
                        -obj.offset_x[0]
                        + (obj.horizontal_movable_size[0] - obj.vertical_movable_size[0]) / 2
                    ),
                    obj.interpiece_offset[1],
                    direction_setting["handle_z_direction"]
                    * (
                    obj.door_size[2] / 2
                        + obj.fixed_part_size[2]
                        + obj.vertical_movable_size[2]
                        + obj.horizontal_movable_size[2] / 2
                    ),
                ])
                top_mesh_rotation = np.array([0, 0, 0])
                __pt = inverse_transformation(__pt, top_mesh_position, top_mesh_rotation)

                if (__pt[0] >= -obj.horizontal_movable_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= obj.horizontal_movable_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.horizontal_movable_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.horizontal_movable_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    __pt[2] >= -obj.horizontal_movable_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= obj.horizontal_movable_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True

            return False
        


        elif (isinstance(obj, PiShape_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), offset_first=True)[0]

            for direction_setting in obj.direction_settings:
                __pt = apply_transformation(
                    np.array([_pt]), 
                    -np.array([direction_setting["handle_y_axis"], 0, 0]),
                    -np.array([0, direction_setting["handle_y_rotation"], 0]),
                    offset_first=True
                )[0]
                __pt = apply_transformation(
                    np.array([__pt]), 
                    -np.array([-direction_setting["handle_y_axis"], 0, 0]),
                    np.array([0, 0, 0]),
                    offset_first=True
                )[0]

                main_mesh_position = np.array([
                    direction_setting["handle_x_direction"] * obj.offset_x[0] + direction_setting["handle_x_position"],
                    obj.interpiece_offset[0] + obj.separation[0] / 2 + obj.sub_size[1] / 2,
                    direction_setting["handle_z_direction"] * (obj.door_size[2] / 2 + obj.sub_size[2] + obj.main_size[2] / 2),
                ])
                main_mesh_rotation = np.array([0, 0, 0])
                __pt = inverse_transformation(__pt, main_mesh_position, main_mesh_rotation)

                if (__pt[0] >= -obj.main_size[0]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[0] <= obj.main_size[0]/2 + AFFORDACE_PROXIMITY_THRES and 
                    __pt[1] >= -obj.main_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.main_size[1]/2 + AFFORDACE_PROXIMITY_THRES and
                    __pt[2] >= -obj.main_size[2]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[2] <= obj.main_size[2]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True

            return False
        

        elif (isinstance(obj, Cylindrical_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), offset_first=True)[0]

            for direction_setting in obj.direction_settings:
                __pt = apply_transformation(
                    np.array([_pt]), 
                    -np.array([direction_setting["handle_y_axis"], 0, 0]),
                    -np.array([0, direction_setting["handle_y_rotation"], 0]),
                    offset_first=True
                )[0]
                __pt = apply_transformation(
                    np.array([__pt]), 
                    -np.array([-direction_setting["handle_y_axis"], 0, 0]),
                    np.array([0, 0, 0]),
                    offset_first=True
                )[0]

                main_mesh_position = np.array([
                    -direction_setting["handle_x_direction"] * obj.offset_x[0] + direction_setting["handle_x_position"],
                    0,
                    direction_setting["handle_z_direction"] * (obj.door_size[2] / 2 + obj.fixed_part_size[1] + obj.sub_size[1] + obj.main_size[1] / 2)
                ])
                main_mesh_rotation = np.array([np.pi/2, 0, 0])
                __pt = inverse_transformation(__pt, main_mesh_position, main_mesh_rotation)

                _pt_radius = np.linalg.norm(__pt[[0, 2]])

                if (_pt_radius <= obj.main_size[0] + AFFORDACE_PROXIMITY_THRES and  
                    __pt[1] >= -obj.main_size[1]/2 - AFFORDACE_PROXIMITY_THRES and
                    __pt[1] <= obj.main_size[1]/2 + AFFORDACE_PROXIMITY_THRES):
                    return True

            return False
        

        elif (isinstance(obj, Spherical_Handle)):
            _pt = apply_transformation(np.array([pt]), -np.array(obj.position), -np.array(obj.rotation), offset_first=True)[0]

            for direction_setting in obj.direction_settings:
                __pt = apply_transformation(
                    np.array([_pt]), 
                    -np.array([direction_setting["handle_y_axis"], 0, 0]),
                    -np.array([0, direction_setting["handle_y_rotation"], 0]),
                    offset_first=True
                )[0]
                __pt = apply_transformation(
                    np.array([__pt]), 
                    -np.array([-direction_setting["handle_y_axis"], 0, 0]),
                    np.array([0, 0, 0]),
                    offset_first=True
                )[0]

                main_mesh_position = np.array([
                    -direction_setting["handle_x_direction"] * obj.offset_x[0] + direction_setting["handle_x_position"],
                    0,
                    direction_setting["handle_z_direction"] * (obj.door_size[2] / 2 + obj.fixed_part_size[1] + obj.sub_size[1] + obj.main_size[1] / 2),
                ])
                main_mesh_rotation = np.array([0, 0, 0])
                __pt = inverse_transformation(__pt, main_mesh_position, main_mesh_rotation)

                _pt_radius = np.linalg.norm(__pt)

                if (_pt_radius <= min(obj.main_size) + AFFORDACE_PROXIMITY_THRES):
                    return True

            return False
    
        
        else:
            return False
    
    
    return is_affordance(obj, pt)


def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT
