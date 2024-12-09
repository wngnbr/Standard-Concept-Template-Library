import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh


class Standard_Door(ConceptTemplate):
    def __init__(self, existence_of_door, size, door_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        door_rotation = [x / 180 * np.pi for x in door_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.existence_of_door = existence_of_door
        self.size = size
        self.door_rotation = door_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        if existence_of_door[0] and existence_of_door[1]:
            left_mesh_position = [size[0] / 2, 0, 0]
            self.left_mesh = Cuboid(size[1], size[0], size[2],
                                    position=left_mesh_position)
            self.left_mesh.vertices = apply_transformation(self.left_mesh.vertices, rotation=[0, -door_rotation[0], 0], position=[-size[0], 0, 0])
            vertices_list.append(self.left_mesh.vertices)
            faces_list.append(self.left_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.left_mesh.vertices)

            self.right_mesh = Cuboid(size[1], size[0], size[2],
                                     position=[-size[0] / 2, 0, 0])
            self.right_mesh.vertices = apply_transformation(self.right_mesh.vertices, rotation=[0, door_rotation[1], 0], position=[size[0], 0, 0])
            vertices_list.append(self.right_mesh.vertices)
            faces_list.append(self.right_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.right_mesh.vertices)

        elif existence_of_door[0]:
            left_mesh_position = [size[0] / 2, 0, 0]
            self.left_mesh = Cuboid(size[1], size[0], size[2],
                                    position=left_mesh_position)
            self.left_mesh.vertices = apply_transformation(self.left_mesh.vertices, rotation=[0, -door_rotation[0], 0], position=[-size[0] / 2, 0, 0])
            vertices_list.append(self.left_mesh.vertices)
            faces_list.append(self.left_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.left_mesh.vertices)

        elif existence_of_door[1]:
            self.right_mesh = Cuboid(size[1], size[0], size[2],
                                     position=[-size[0] / 2, 0, 0])
            self.right_mesh.vertices = apply_transformation(self.right_mesh.vertices, rotation=[0, door_rotation[1], 0], position=[size[0] / 2, 0, 0])
            vertices_list.append(self.right_mesh.vertices)
            faces_list.append(self.right_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.right_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Door'


class Standard_Doorframe(ConceptTemplate):
    def __init__(self, door_size, existence_of_doorframe, main_outer_size, main_inner_outer_offset, main_offset, sub1_outer_size, sub1_inner_size, sub1_inner_outer_offset, sub1_offset,sub2_outer_size, sub2_inner_size, sub2_inner_outer_offset, sub2_offset, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.door_size = door_size
        self.existence_of_doorframe = existence_of_doorframe
        self.main_outer_size = main_outer_size
        self.main_inner_outer_offset = main_inner_outer_offset
        self.main_offset = main_offset
        self.sub1_outer_size = sub1_outer_size
        self.sub1_inner_size = sub1_inner_size
        self.sub1_inner_outer_offset = sub1_inner_outer_offset
        self.sub1_offset = sub1_offset
        self.sub2_outer_size = sub2_outer_size
        self.sub2_inner_size = sub2_inner_size
        self.sub2_inner_outer_offset = sub2_inner_outer_offset
        self.sub2_offset = sub2_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        main_top_size = [
            door_size[0],
            main_outer_size[1] + main_inner_outer_offset[1] - door_size[1],
            main_outer_size[2],
        ]
        main_left_size = [
            (main_outer_size[0] - door_size[0]) / 2 - main_inner_outer_offset[0],
            main_outer_size[1],
            main_outer_size[2],
        ]
        main_right_size = [
            (main_outer_size[0] - door_size[0]) / 2 + main_inner_outer_offset[0],
            main_outer_size[1],
            main_outer_size[2],
        ]

        main_top_mesh_position = [
            main_offset[0],
            door_size[1] / 2 + main_top_size[1] / 2,
            main_offset[1],
        ]
        self.main_top_mesh = Cuboid(main_top_size[1], main_top_size[0], main_top_size[2],
                                    position=main_top_mesh_position)
        vertices_list.append(self.main_top_mesh.vertices)
        faces_list.append(self.main_top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.main_top_mesh.vertices)

        main_left_mesh_position = [
            main_offset[0] - door_size[0] / 2 - main_left_size[0] / 2,
            -door_size[1] / 2 + main_outer_size[1] / 2 + main_inner_outer_offset[1],
            main_offset[1],
        ]
        self.main_left_mesh = Cuboid(main_left_size[1], main_left_size[0], main_left_size[2],
                                     position=main_left_mesh_position)
        vertices_list.append(self.main_left_mesh.vertices)
        faces_list.append(self.main_left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.main_left_mesh.vertices)

        main_right_mesh_position = [
            main_offset[0] + door_size[0] / 2 + main_right_size[0] / 2,
            -door_size[1] / 2 + main_outer_size[1] / 2 + main_inner_outer_offset[1],
            main_offset[1],
        ]
        self.main_right_mesh = Cuboid(main_right_size[1], main_right_size[0], main_right_size[2],
                                      position=main_right_mesh_position)
        vertices_list.append(self.main_right_mesh.vertices)
        faces_list.append(self.main_right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.main_right_mesh.vertices)

        if existence_of_doorframe[0]:
            front_top_size = [
                sub1_inner_size[0],
                sub1_outer_size[1] + sub1_inner_outer_offset[1] - sub1_inner_size[1],
                sub1_outer_size[2],
            ]
            front_left_size = [
                (sub1_outer_size[0] - sub1_inner_size[0]) / 2
                - sub1_inner_outer_offset[0],
                sub1_outer_size[1],
                sub1_outer_size[2],
            ]
            front_right_size = [
                (sub1_outer_size[0] - sub1_inner_size[0]) / 2
                + sub1_inner_outer_offset[0],
                sub1_outer_size[1],
                sub1_outer_size[2],
            ]

            front_top_mesh_position = [
                sub1_offset[0],
                sub1_offset[1] - door_size[1] / 2 + sub1_inner_size[1] + front_top_size[1] / 2,
                -sub1_outer_size[2] / 2 + main_offset[1] + main_outer_size[2] / 2,
            ]
            self.front_top_mesh = Cuboid(front_top_size[1], front_top_size[0], front_top_size[2],
                                         position=front_top_mesh_position)
            vertices_list.append(self.front_top_mesh.vertices)
            faces_list.append(self.front_top_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.front_top_mesh.vertices)

            front_left_mesh_position = [
                sub1_offset[0] - sub1_inner_size[0] / 2 - front_left_size[0] / 2,
                sub1_offset[1] - door_size[1] / 2 + sub1_outer_size[1] / 2 + sub1_inner_outer_offset[1],
                -sub1_outer_size[2] / 2 + main_offset[1] + main_outer_size[2] / 2,
            ]
            self.front_left_mesh = Cuboid(front_left_size[1], front_left_size[0], front_left_size[2],
                                          position=front_left_mesh_position)
            vertices_list.append(self.front_left_mesh.vertices)
            faces_list.append(self.front_left_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.front_left_mesh.vertices)

            front_right_mesh_position = [
                sub1_offset[0] + sub1_inner_size[0] / 2 + front_left_size[0] / 2,
                sub1_offset[1] - door_size[1] / 2 + sub1_outer_size[1] / 2 + sub1_inner_outer_offset[1],
                -sub1_outer_size[2] / 2 + main_offset[1] + main_outer_size[2] / 2,
            ]
            self.front_right_mesh = Cuboid(front_right_size[1], front_right_size[0], front_right_size[2],
                                           position=front_right_mesh_position)
            vertices_list.append(self.front_right_mesh.vertices)
            faces_list.append(self.front_right_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.front_right_mesh.vertices)

        if existence_of_doorframe[1]:
            back_top_size = [
                sub2_inner_size[0],
                sub2_outer_size[1] + sub2_inner_outer_offset[1] - sub2_inner_size[1],
                sub2_outer_size[2],
            ]
            back_left_size = [
                (sub2_outer_size[0] - sub2_inner_size[0]) / 2
                - sub2_inner_outer_offset[0],
                sub2_outer_size[1],
                sub2_outer_size[2],
            ]
            back_right_size = [
                (sub2_outer_size[0] - sub2_inner_size[0]) / 2
                + sub2_inner_outer_offset[0],
                sub2_outer_size[1],
                sub2_outer_size[2],
            ]

            back_top_mesh_position = [
                sub2_offset[0],
                sub2_offset[1] - door_size[1] / 2 + sub2_inner_size[1] + back_top_size[1] / 2,
                -sub2_outer_size[2] / 2 + main_offset[1] - main_outer_size[2] / 2,
            ]
            self.back_top_mesh = Cuboid(back_top_size[1], back_top_size[0], back_top_size[2],
                                        position=back_top_mesh_position)
            vertices_list.append(self.back_top_mesh.vertices)
            faces_list.append(self.back_top_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.back_top_mesh.vertices)

            back_left_mesh_position = [
                sub2_offset[0] - sub2_inner_size[0] / 2 - back_left_size[0] / 2,
                sub2_offset[1] - door_size[1] / 2 + sub2_outer_size[1] / 2 + sub2_inner_outer_offset[1],
                -sub2_outer_size[2] / 2 + main_offset[1] - main_outer_size[2] / 2,
            ]
            self.back_left_mesh = Cuboid(back_left_size[1], back_left_size[0], back_left_size[2],
                                         position=back_left_mesh_position)
            vertices_list.append(self.back_left_mesh.vertices)
            faces_list.append(self.back_left_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.back_left_mesh.vertices)

            back_right_mesh_position = [
                sub2_offset[0] + sub2_inner_size[0] / 2 + back_right_size[0] / 2,
                sub2_offset[1] - door_size[1] / 2 + sub2_outer_size[1] / 2 + sub2_inner_outer_offset[1],
                -sub2_outer_size[2] / 2 + main_offset[1] - main_outer_size[2] / 2,
            ]
            self.back_right_mesh = Cuboid(back_right_size[1], back_right_size[0], back_right_size[2],
                                          position=back_right_mesh_position)
            vertices_list.append(self.back_right_mesh.vertices)
            faces_list.append(self.back_right_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.back_right_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Doorframe'


class Standard_Hinge(ConceptTemplate):
    def __init__(self, existence_of_door, number_of_hinge, size, separation, offset_1, offset_2, position=[0, 0, 0], rotation=[0, 0, 0]):
        
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.existence_of_door = existence_of_door
        self.number_of_hinge = number_of_hinge
        self.size = size
        self.separation = separation
        self.offset_1 = offset_1
        self.offset_2 = offset_2

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        offset = offset_1
        for existence in range(existence_of_door[0] + existence_of_door[1]):
            for i in range(number_of_hinge[0]):
                tmp_mesh_position = [
                    offset[0],
                    offset[1] + size[1] * (i + 0.5) + sum(separation[0:i]),
                    offset[2],
                ]
                self.tmp_mesh = Cylinder(size[1], size[0], size[0],
                                    position=tmp_mesh_position)

                vertices_list.append(self.tmp_mesh.vertices)
                faces_list.append(self.tmp_mesh.faces + total_num_vertices)
                total_num_vertices += len(self.tmp_mesh.vertices)
            
            offset = offset_2

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Hinge'


class LShape_Handle(ConceptTemplate):
    def __init__(self, existence_of_door, door_rotation, door_size, existence_of_handle, fixed_part_size, vertical_movable_size, horizontal_movable_size, interpiece_offset, offset_x, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        door_rotation = [x / 180 * np.pi for x in door_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.existence_of_handle = existence_of_handle
        self.door_size = door_size
        self.existence_of_door = existence_of_door
        self.door_rotation = door_rotation
        self.fixed_part_size = fixed_part_size
        self.vertical_movable_size = vertical_movable_size
        self.horizontal_movable_size = horizontal_movable_size
        self.interpiece_offset = interpiece_offset
        self.offset_x = offset_x

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.direction_settings = []

        double_door = 0
        if existence_of_door[0] and existence_of_door[1]:
            double_door = 1

        # parameter calculate
        for door in [0, 1]:
            if not existence_of_door[door]:
                continue
            for handle in [0, 1]:
                if not existence_of_handle[handle]:
                    continue
                # right_door(x+)
                if door:
                    handle_y_rotation = door_rotation[1]
                    handle_x_direction = 1
                    handle_y_axis = door_size[0] / 2
                # left_door(x-)
                else:
                    handle_y_rotation = -door_rotation[0]
                    handle_x_direction = -1
                    handle_y_axis = -door_size[0] / 2
                # front_handle(z+)
                if handle:
                    handle_z_direction = 1
                # back_handle(z-)
                else:
                    handle_z_direction = -1
                if double_door:
                    handle_x_position = handle_y_axis
                    handle_y_axis *= 2
                else:
                    handle_x_position = 0
                self.direction_settings.append(
                    {
                        "handle_x_direction": handle_x_direction,
                        "handle_x_position": handle_x_position,
                        "handle_z_direction": handle_z_direction,
                        "handle_y_axis": handle_y_axis,
                        "handle_y_rotation": handle_y_rotation,
                    }
                )

        # meshes definition
        for direction_setting in self.direction_settings:
            tmp_meshes = []

            base_mesh_position = [
                -direction_setting["handle_x_direction"] * offset_x[0] + direction_setting["handle_x_position"],
                0,
                direction_setting["handle_z_direction"] * (door_size[2] / 2 + fixed_part_size[2] / 2),
            ]
            self.base_mesh = Cuboid(fixed_part_size[1], fixed_part_size[0], fixed_part_size[2],
                               position=base_mesh_position)
            tmp_meshes.append(self.base_mesh)

            middle_mesh_position = [
                direction_setting["handle_x_position"] + interpiece_offset[0] - direction_setting["handle_x_direction"] * offset_x[0],
                interpiece_offset[1],
                direction_setting["handle_z_direction"] * (door_size[2] / 2 + fixed_part_size[2] + vertical_movable_size[2] / 2)
            ]
            self.middle_mesh = Cuboid(vertical_movable_size[1], vertical_movable_size[0], vertical_movable_size[2],
                                 position=middle_mesh_position)
            tmp_meshes.append(self.middle_mesh)

            top_mesh_position = [
                direction_setting["handle_x_position"]
                + interpiece_offset[0]
                + direction_setting["handle_x_direction"]
                * (
                    -offset_x[0]
                    + (horizontal_movable_size[0] - vertical_movable_size[0]) / 2
                ),
                interpiece_offset[1],
                direction_setting["handle_z_direction"]
                * (
                door_size[2] / 2
                    + fixed_part_size[2]
                    + vertical_movable_size[2]
                    + horizontal_movable_size[2] / 2
                ),
            ]
            self.top_mesh = Cuboid(horizontal_movable_size[1], horizontal_movable_size[0], horizontal_movable_size[2],
                              position=top_mesh_position)
            tmp_meshes.append(self.top_mesh)

            for tmp_mesh in tmp_meshes:
                tmp_mesh.vertices = apply_transformation(tmp_mesh.vertices, [-direction_setting["handle_y_axis"], 0, 0], [0, 0, 0])
                tmp_mesh.vertices = apply_transformation(tmp_mesh.vertices, [direction_setting["handle_y_axis"], 0, 0], [0, direction_setting["handle_y_rotation"], 0])
                vertices_list.append(tmp_mesh.vertices)
                faces_list.append(tmp_mesh.faces + total_num_vertices)
                total_num_vertices += len(tmp_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class PiShape_Handle(ConceptTemplate):
    def __init__(self, existence_of_door, door_rotation, door_size, existence_of_handle, main_size, sub_size, separation, interpiece_offset, offset_x, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        door_rotation = [x / 180 * np.pi for x in door_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.existence_of_handle = existence_of_handle
        self.door_size = door_size
        self.existence_of_door = existence_of_door
        self.door_rotation = door_rotation
        self.main_size = main_size
        self.sub_size = sub_size
        self.separation = separation
        self.interpiece_offset = interpiece_offset
        self.offset_x = offset_x

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.direction_settings = []
        double_door = 0
        if existence_of_door[0] and existence_of_door[1]:
            double_door = 1

        # parameter calculate
        for door in [0, 1]:
            if not existence_of_door[door]:
                continue
            for handle in [0, 1]:
                if not existence_of_handle[handle]:
                    continue
                # right_door(x+)
                if door:
                    handle_y_rotation = door_rotation[1]
                    handle_x_direction = 1
                    handle_y_axis = door_size[0] / 2
                # left_door(x-)
                else:
                    handle_y_rotation = -door_rotation[0]
                    handle_x_direction = -1
                    handle_y_axis = -door_size[0] / 2
                # front_handle(z+)
                if handle:
                    handle_z_direction = -1
                # back_handle(z-)
                else:
                    handle_z_direction = 1
                if double_door:
                    handle_x_position = handle_y_axis
                    handle_y_axis *= 2
                else:
                    handle_x_position = 0
                self.direction_settings.append(
                    {
                        "handle_x_direction": handle_x_direction,
                        "handle_x_position": handle_x_position,
                        "handle_z_direction": handle_z_direction,
                        "handle_y_axis": handle_y_axis,
                        "handle_y_rotation": handle_y_rotation,
                    }
                )

        # meshes definition
        for direction_setting in self.direction_settings:
            tmp_meshes = []

            top_mesh_position = [
                direction_setting["handle_x_direction"] * offset_x[0] + direction_setting["handle_x_position"],
                separation[0] + sub_size[1],
                direction_setting["handle_z_direction"] * (door_size[2] / 2 + sub_size[2] / 2),
            ]
            self.top_mesh = Cuboid(sub_size[1], sub_size[0], sub_size[2],
                              position=top_mesh_position)
            tmp_meshes.append(self.top_mesh)

            bottom_mesh_position = [
                direction_setting["handle_x_direction"] * offset_x[0] + direction_setting["handle_x_position"],
                0,
                direction_setting["handle_z_direction"] * (door_size[2] / 2 + sub_size[2] / 2),
            ]
            self.bottom_mesh = Cuboid(sub_size[1], sub_size[0], sub_size[2],
                                 position=bottom_mesh_position)
            tmp_meshes.append(self.bottom_mesh)

            main_mesh_position = [
                direction_setting["handle_x_direction"] * offset_x[0] + direction_setting["handle_x_position"],
                interpiece_offset[0] + separation[0] / 2 + sub_size[1] / 2,
                direction_setting["handle_z_direction"] * (door_size[2] / 2 + sub_size[2] + main_size[2] / 2),
            ]
            self.main_mesh = Cuboid(main_size[1], main_size[0], main_size[2],
                               position=main_mesh_position)
            tmp_meshes.append(self.main_mesh)

            for tmp_mesh in tmp_meshes:
                tmp_mesh.vertices = apply_transformation(tmp_mesh.vertices, [-direction_setting["handle_y_axis"], 0, 0], [0, 0, 0])
                tmp_mesh.vertices = apply_transformation(tmp_mesh.vertices, [direction_setting["handle_y_axis"], 0, 0], [0, direction_setting["handle_y_rotation"], 0])
                vertices_list.append(tmp_mesh.vertices)
                faces_list.append(tmp_mesh.faces + total_num_vertices)
                total_num_vertices += len(tmp_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Cylindrical_Handle(ConceptTemplate):
    def __init__(self, existence_of_door, door_rotation, door_size, existence_of_handle, fixed_part_size, sub_size, main_size, offset_x, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        door_rotation = [x / 180 * np.pi for x in door_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.existence_of_handle = existence_of_handle
        self.door_size = door_size
        self.existence_of_door = existence_of_door
        self.door_rotation = door_rotation
        self.main_size = main_size
        self.sub_size = sub_size
        self.fixed_part_size = fixed_part_size
        self.offset_x = offset_x

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.direction_settings = []
        double_door = 0
        if existence_of_door[0] and existence_of_door[1]:
            double_door = 1

        # parameter calculate
        for door in [0, 1]:
            if not existence_of_door[door]:
                continue
            for handle in [0, 1]:
                if not existence_of_handle[handle]:
                    continue
                # right_door(x+)
                if door:
                    handle_y_rotation = door_rotation[1]
                    handle_x_direction = 1
                    handle_y_axis = door_size[0] / 2
                # left_door(x-)
                else:
                    handle_y_rotation = -door_rotation[0]
                    handle_x_direction = -1
                    handle_y_axis = -door_size[0] / 2
                # front_handle(z+)
                if handle:
                    handle_z_direction = 1
                # back_handle(z-)
                else:
                    handle_z_direction = -1
                if double_door:
                    handle_x_position = handle_y_axis
                    handle_y_axis *= 2
                else:
                    handle_x_position = 0
                self.direction_settings.append(
                    {
                        "handle_x_direction": handle_x_direction,
                        "handle_x_position": handle_x_position,
                        "handle_z_direction": handle_z_direction,
                        "handle_y_axis": handle_y_axis,
                        "handle_y_rotation": handle_y_rotation,
                    }
                )

        # meshes definition
        for direction_setting in self.direction_settings:
            tmp_meshes = []

            base_mesh_position = [
                -direction_setting["handle_x_direction"] * offset_x[0] + direction_setting["handle_x_position"],
                0,
                direction_setting["handle_z_direction"] * (door_size[2] / 2 + fixed_part_size[1] / 2),
            ]
            base_mesh_rotation = [np.pi / 2, 0, 0]
            self.base_mesh = Cylinder(fixed_part_size[1], fixed_part_size[0], fixed_part_size[0],
                                 position=base_mesh_position,
                                 rotation=base_mesh_rotation)
            tmp_meshes.append(self.base_mesh)

            middle_mesh_position = [
                -direction_setting["handle_x_direction"] * offset_x[0] + direction_setting["handle_x_position"],
                0,
                direction_setting["handle_z_direction"] * (door_size[2] / 2 + fixed_part_size[1] + sub_size[1] / 2),
            ]
            middle_mesh_rotation = [np.pi / 2, 0, 0]
            self.middle_mesh = Cylinder(sub_size[1], sub_size[0], sub_size[0],
                                   position=middle_mesh_position,
                                   rotation=middle_mesh_rotation)
            tmp_meshes.append(self.middle_mesh)

            main_mesh_position = [
                -direction_setting["handle_x_direction"] * offset_x[0] + direction_setting["handle_x_position"],
                0,
                direction_setting["handle_z_direction"] * (door_size[2] / 2 + fixed_part_size[1] + sub_size[1] + main_size[1] / 2)
            ]
            main_mesh_rotation = [np.pi / 2, 0, 0]
            self.main_mesh = Cylinder(main_size[1], main_size[0], main_size[0],
                                 position=main_mesh_position,
                                 rotation=main_mesh_rotation)
            tmp_meshes.append(self.main_mesh)

            for tmp_mesh in tmp_meshes:
                tmp_mesh.vertices = apply_transformation(tmp_mesh.vertices, [-direction_setting["handle_y_axis"], 0, 0], [0, 0, 0])
                tmp_mesh.vertices = apply_transformation(tmp_mesh.vertices, [direction_setting["handle_y_axis"], 0, 0], [0, direction_setting["handle_y_rotation"], 0])
                vertices_list.append(tmp_mesh.vertices)
                faces_list.append(tmp_mesh.faces + total_num_vertices)
                total_num_vertices += len(tmp_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Spherical_Handle(ConceptTemplate):
    def __init__(self, existence_of_door, door_rotation, door_size, existence_of_handle, fixed_part_size, sub_size, main_size, offset_x, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        door_rotation = [x / 180 * np.pi for x in door_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.existence_of_handle = existence_of_handle
        self.door_size = door_size
        self.existence_of_door = existence_of_door
        self.door_rotation = door_rotation
        self.main_size = main_size
        self.sub_size = sub_size
        self.fixed_part_size = fixed_part_size
        self.offset_x = offset_x

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.direction_settings = []
        double_door = 0
        if existence_of_door[0] and existence_of_door[1]:
            double_door = 1

        # parameter calculate
        for door in [0, 1]:
            if not existence_of_door[door]:
                continue
            for handle in [0, 1]:
                if not existence_of_handle[handle]:
                    continue
                # right_door(x+)
                if door:
                    handle_y_rotation = door_rotation[1]
                    handle_x_direction = 1
                    handle_y_axis = door_size[0] / 2
                # left_door(x-)
                else:
                    handle_y_rotation = -door_rotation[0]
                    handle_x_direction = -1
                    handle_y_axis = -door_size[0] / 2
                # front_handle(z+)
                if handle:
                    handle_z_direction = 1
                # back_handle(z-)
                else:
                    handle_z_direction = -1
                if double_door:
                    handle_x_position = handle_y_axis
                    handle_y_axis *= 2
                else:
                    handle_x_position = 0
                self.direction_settings.append(
                    {
                        "handle_x_direction": handle_x_direction,
                        "handle_x_position": handle_x_position,
                        "handle_z_direction": handle_z_direction,
                        "handle_y_axis": handle_y_axis,
                        "handle_y_rotation": handle_y_rotation,
                    }
                )
        # meshes definition
        for direction_setting in self.direction_settings:
            tmp_meshes = []

            base_mesh_position = [
                -direction_setting["handle_x_direction"] * offset_x[0] + direction_setting["handle_x_position"],
                0,
                direction_setting["handle_z_direction"] * (door_size[2] / 2 + fixed_part_size[1] / 2),
            ]
            base_mesh_rotation = [np.pi / 2, 0, 0]
            self.base_mesh = Cylinder(fixed_part_size[1], fixed_part_size[0], fixed_part_size[0],
                                      position=base_mesh_position,
                                      rotation=base_mesh_rotation)
            tmp_meshes.append(self.base_mesh)

            middle_mesh_position = [
                -direction_setting["handle_x_direction"] * offset_x[0] + direction_setting["handle_x_position"],
                0,
                direction_setting["handle_z_direction"] * (door_size[2] / 2 + fixed_part_size[1] + sub_size[1] / 2),
            ]
            middle_mesh_rotation = [np.pi / 2, 0, 0]
            self.middle_mesh = Cylinder(sub_size[1], sub_size[0], sub_size[0],
                                        position=middle_mesh_position,
                                        rotation=middle_mesh_rotation)
            tmp_meshes.append(self.middle_mesh)

            main_mesh_position = [
                -direction_setting["handle_x_direction"] * offset_x[0] + direction_setting["handle_x_position"],
                0,
                direction_setting["handle_z_direction"] * (door_size[2] / 2 + fixed_part_size[1] + sub_size[1] + main_size[1] / 2),
            ]
            self.main_mesh = Sphere(main_size[0], radius_y=main_size[1], radius_z=main_size[2],
                                    position=main_mesh_position)
            tmp_meshes.append(self.main_mesh)

            for tmp_mesh in tmp_meshes:
                tmp_mesh.vertices = apply_transformation(tmp_mesh.vertices, [-direction_setting["handle_y_axis"], 0, 0], [0, 0, 0])
                tmp_mesh.vertices = apply_transformation(tmp_mesh.vertices, [direction_setting["handle_y_axis"], 0, 0], [0, direction_setting["handle_y_rotation"], 0])
                vertices_list.append(tmp_mesh.vertices)
                faces_list.append(tmp_mesh.faces + total_num_vertices)
                total_num_vertices += len(tmp_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)
        
        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'