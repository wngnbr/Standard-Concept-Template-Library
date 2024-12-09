import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation, get_rodrigues_matrix
from knowledge_utils import *
import trimesh


class Regular_seat(ConceptTemplate):
    def __init__(self, size, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.seat_mesh = Cuboid(self.size[1], self.size[0], self.size[2])
        vertices_list.append(self.seat_mesh.vertices)
        faces_list.append(self.seat_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.seat_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Seat'


class Round_seat(ConceptTemplate):
    def __init__(self, size, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.radius = size[0]
        self.height = size[1]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.seat_mesh = Cylinder(self.height, self.radius)
        vertices_list.append(self.seat_mesh.vertices)
        faces_list.append(self.seat_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.seat_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Seat'


class Solid_back(ConceptTemplate):
    def __init__(self, size, back_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        back_rotation = [x / 180 * np.pi for x in back_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.back_rotation = back_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_rotation = [back_rotation[0], 0, 0]
        mesh_position = [0, size[1] / 2 * np.cos(back_rotation[0]), 0]
        self.back_mesh = Cuboid(self.size[1], self.size[0], self.size[2], 
                                position=mesh_position, rotation=mesh_rotation)
        vertices_list.append(self.back_mesh.vertices)
        faces_list.append(self.back_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.back_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Back'


class Ladder_back(ConceptTemplate):
    def __init__(self, main_horizontal_piece_size, main_vertical_piece_size, sub_horizontal_piece_size,
                 main_vertical_separation, sub_offset, interval_between_subs, back_rotation, number_of_subs,
                 position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        back_rotation = [x / 180 * np.pi for x in back_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.main_horizontal_piece_size = main_horizontal_piece_size
        self.main_vertical_piece_size = main_vertical_piece_size
        self.sub_horizontal_piece_size = sub_horizontal_piece_size
        self.main_vertical_separation = main_vertical_separation
        self.sub_offset = sub_offset
        self.interval_between_subs = interval_between_subs
        self.back_rotation = back_rotation
        self.number_of_subs = number_of_subs

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(3 + number_of_subs[0]):
            mesh_rotation = [back_rotation[0], 0, 0]
            if i < 2:
                flag = 1 if i == 1 else -1
                vertical_position = [flag * main_vertical_separation[0] / 2,
                                     main_vertical_piece_size[1] / 2 * np.cos(back_rotation[0]), 0]
                self.mesh = Cuboid(self.main_vertical_piece_size[1], self.main_vertical_piece_size[0],
                                   self.main_vertical_piece_size[2],
                                   position=vertical_position, rotation=mesh_rotation)
            elif i == 2:
                horizontal_position = [0,
                                       (main_vertical_piece_size[1] + main_horizontal_piece_size[1] / 2) * np.cos(
                                           back_rotation[0]),
                                       (main_vertical_piece_size[1] + main_horizontal_piece_size[1]) / 2 * np.sin(
                                           back_rotation[0])]
                self.mesh = Cuboid(self.main_horizontal_piece_size[1], self.main_horizontal_piece_size[0],
                                   self.main_horizontal_piece_size[2],
                                   position=horizontal_position, rotation=mesh_rotation)
            else:
                sub_position = [0, (main_vertical_piece_size[1] / 2 + sub_offset[0] - (i - 3) * interval_between_subs[
                    0]) * np.cos(back_rotation[0]),
                                (sub_offset[0] - (i - 3) * interval_between_subs[0]) * np.sin(back_rotation[0])]
                self.mesh = Cuboid(self.sub_horizontal_piece_size[0],
                                   self.main_vertical_separation[0] - self.main_vertical_piece_size[0],
                                   self.sub_horizontal_piece_size[1],
                                   position=sub_position, rotation=mesh_rotation)

            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Back' 


class Splat_back(ConceptTemplate):
    def __init__(self, main_horizontal_piece_size, main_vertical_piece_size, sub_vertical_piece_size,
                 main_vertical_separation, sub_offset, interval_between_subs, back_rotation, number_of_subs,
                 position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        back_rotation = [x / 180 * np.pi for x in back_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.main_horizontal_piece_size = main_horizontal_piece_size
        self.main_vertical_piece_size = main_vertical_piece_size
        self.sub_vertical_piece_size = sub_vertical_piece_size
        self.main_vertical_separation = main_vertical_separation
        self.sub_offset = sub_offset
        self.interval_between_subs = interval_between_subs
        self.back_rotation = back_rotation
        self.number_of_subs = number_of_subs

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(3 + number_of_subs[0]):
            mesh_rotation = [back_rotation[0], 0, 0]
            if i < 2:
                flag = 1 if i == 1 else -1
                vertical_position = [flag * main_vertical_separation[0] / 2,
                                     main_vertical_piece_size[1] / 2 * np.cos(back_rotation[0]), 0]
                self.mesh = Cuboid(self.main_vertical_piece_size[1], self.main_vertical_piece_size[0],
                                   self.main_vertical_piece_size[2],
                                   position=vertical_position, rotation=mesh_rotation)
            elif i == 2:
                horizontal_position = [0,
                                       (main_vertical_piece_size[1] + main_horizontal_piece_size[1] / 2) * np.cos(
                                           back_rotation[0]),
                                       (main_vertical_piece_size[1] + main_horizontal_piece_size[1]) / 2 * np.sin(
                                           back_rotation[0])]
                self.mesh = Cuboid(self.main_horizontal_piece_size[1], self.main_horizontal_piece_size[0],
                                   self.main_horizontal_piece_size[2],
                                   position=horizontal_position, rotation=mesh_rotation)
            else:
                sub_position = [sub_offset[0] + (i - 3) * interval_between_subs[0],
                                main_vertical_piece_size[1] / 2 * np.cos(back_rotation[0]), 0]
                self.mesh = Cuboid(self.main_vertical_piece_size[1],
                                   self.sub_vertical_piece_size[0],
                                   self.sub_vertical_piece_size[1],
                                   position=sub_position, rotation=mesh_rotation)

            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Back'


class Latice_back(ConceptTemplate):
    def __init__(self, main_horizontal_piece_size, main_vertical_piece_size, main_vertical_separation,
                 sub_vertical_piece_size, sub_horizontal_piece_size, sub_horizontal_offset, sub_vertical_offset,
                 interval_between_subs, back_rotation, number_of_subs, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        back_rotation = [x / 180 * np.pi for x in back_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.main_horizontal_piece_size = main_horizontal_piece_size
        self.main_vertical_piece_size = main_vertical_piece_size
        self.sub_vertical_piece_size = sub_vertical_piece_size
        self.sub_horizontal_piece_size = sub_horizontal_piece_size
        self.main_vertical_separation = main_vertical_separation
        self.sub_horizontal_offset = sub_horizontal_offset
        self.sub_vertical_offset = sub_vertical_offset
        self.interval_between_subs = interval_between_subs
        self.back_rotation = back_rotation
        self.number_of_subs = number_of_subs

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(4 + number_of_subs[0]):
            mesh_rotation = [back_rotation[0], 0, 0]
            if i < 2:
                flag = 1 if i == 1 else -1
                vertical_position = [flag * main_vertical_separation[0] / 2,
                                     main_vertical_piece_size[1] / 2 * np.cos(back_rotation[0]), 0]
                self.mesh = Cuboid(self.main_vertical_piece_size[1], self.main_vertical_piece_size[0],
                                   self.main_vertical_piece_size[2],
                                   position=vertical_position, rotation=mesh_rotation)
            elif i == 2:
                horizontal_position = [0,
                                       (main_vertical_piece_size[1] + main_horizontal_piece_size[1] / 2) * np.cos(
                                           back_rotation[0]),
                                       (main_vertical_piece_size[1] + main_horizontal_piece_size[1]) / 2 * np.sin(
                                           back_rotation[0])]
                self.mesh = Cuboid(self.main_horizontal_piece_size[1], self.main_horizontal_piece_size[0],
                                   self.main_horizontal_piece_size[2],
                                   position=horizontal_position, rotation=mesh_rotation)
            elif i == 3:
                sub_horizontal_position = [0,
                                           (main_vertical_piece_size[1] / 2 + sub_horizontal_offset[0]) * np.cos(
                                               back_rotation[0]),
                                           sub_horizontal_offset[0] * np.sin(back_rotation[0])]
                self.mesh = Cuboid(self.sub_horizontal_piece_size[0],
                                   self.main_vertical_separation[0] - self.main_vertical_piece_size[0],
                                   self.sub_horizontal_piece_size[1],
                                   position=sub_horizontal_position, rotation=mesh_rotation)
            else:
                sub_vertical_position = [sub_vertical_offset[0] + (i - 4) * interval_between_subs[0],
                                         (main_vertical_piece_size[1] * 3 / 2 + sub_horizontal_offset[0] +
                                          sub_horizontal_piece_size[0] / 2) / 2 * np.cos(back_rotation[0]),
                                         (main_vertical_piece_size[1] / 2 + sub_horizontal_offset[0] +
                                          sub_horizontal_piece_size[0] / 2) / 2 * np.sin(back_rotation[0])]
                self.mesh = Cuboid(
                    self.main_vertical_piece_size[1] / 2 - sub_horizontal_piece_size[0] / 2 - sub_horizontal_offset[0],
                    self.sub_vertical_piece_size[0],
                    self.sub_vertical_piece_size[1],
                    position=sub_vertical_position, rotation=mesh_rotation)

            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)


        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Back'


class Slat_back(ConceptTemplate):
    def __init__(self, main_vertical_piece_size, sub_horizontal_piece_size, main_vertical_separation,
                 sub_horizontal_offset, interval_between_subs, main_vertical_rotation, back_rotation, number_of_subs,
                 position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        back_rotation = [x / 180 * np.pi for x in back_rotation]
        main_vertical_rotation = [x / 180 * np.pi for x in main_vertical_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.main_vertical_piece_size = main_vertical_piece_size
        self.sub_horizontal_piece_size = sub_horizontal_piece_size
        self.main_vertical_separation = main_vertical_separation
        self.sub_horizontal_offset = sub_horizontal_offset
        self.interval_between_subs = interval_between_subs
        self.main_vertical_rotation = main_vertical_rotation
        self.back_rotation = back_rotation
        self.number_of_subs = number_of_subs

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(2 + number_of_subs[0]):
            if i < 2:
                flag = 1 if i == 1 else -1
                vertical_rotation = [back_rotation[0], flag * main_vertical_rotation[0], 0]
                vertical_position = [flag * main_vertical_separation[0] / 2,
                                     main_vertical_piece_size[1] / 2 * np.cos(back_rotation[0]), 0]
                self.mesh = Cuboid(self.main_vertical_piece_size[1], self.main_vertical_piece_size[0],
                                   self.main_vertical_piece_size[2],
                                   position=vertical_position, rotation=vertical_rotation, rotation_order="ZYX")
            else:
                sub_horizontal_rotation = [back_rotation[0], 0, 0]
                sub_horizontal_position = [0,
                                           (main_vertical_piece_size[1] / 2 + sub_horizontal_offset[0] -
                                            (i - 2) * interval_between_subs[0]) * np.cos(back_rotation[0]) -
                                           main_vertical_piece_size[0] / 2 * np.sin(main_vertical_rotation[0]) *
                                           np.sin(back_rotation[0]),
                                           (sub_horizontal_offset[0] - (i - 2) * interval_between_subs[0]) *
                                           np.sin(back_rotation[0]) + main_vertical_piece_size[0] / 2 *
                                           np.sin(main_vertical_rotation[0]) * np.cos(back_rotation[0])]
                self.mesh = Cuboid(
                    self.sub_horizontal_piece_size[0],
                    self.main_vertical_separation[0] - main_vertical_piece_size[0] * np.cos(main_vertical_rotation[0]),
                    self.sub_horizontal_piece_size[1],
                    position=sub_horizontal_position, rotation=sub_horizontal_rotation)

            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Back'


class Regular_leg(ConceptTemplate):
    def __init__(self, number_of_legs, legs_separation, central_rotation,
                 symmetry_mode, front_legs_size, front_rotation, rear_legs_size, rear_rotation,
                 position=[0, 0, 0], rotation=[0, 0, 0]):
        
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        central_rotation = [x / 180 * np.pi for x in central_rotation]
        front_rotation = [x / 180 * np.pi for x in front_rotation]
        rear_rotation = [x / 180 * np.pi for x in rear_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_legs = number_of_legs
        self.legs_separation = legs_separation
        self.central_rotation = central_rotation
        self.symmetry_mode = symmetry_mode
        self.front_legs_size = front_legs_size
        self.front_rotation = front_rotation
        self.rear_legs_size = rear_legs_size
        self.rear_rotation = rear_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(number_of_legs[0]):
            rotation_sign = 1 if (i % 2 == 0 or self.symmetry_mode[0] == 0) else -1
            position_sign = 1 if (i % 2 == 1) else -1
            if number_of_legs[0] == 1:
                mesh_rotation = [front_rotation[0], central_rotation[0], 0]
                mesh_position = [0, -front_legs_size[1] / 2 * np.cos(front_rotation[0]), 0]
                self.mesh = Cuboid(self.front_legs_size[1], self.front_legs_size[0], self.front_legs_size[2],
                                   position=mesh_position, rotation=mesh_rotation)
            if number_of_legs[0] == 2:
                mesh_rotation = [front_rotation[0], central_rotation[0], rotation_sign * front_rotation[1]]
                mesh_position = [position_sign * legs_separation[0] / 2 * np.cos(central_rotation[0]),
                                 -front_legs_size[1] / 2 * np.cos(front_rotation[0]) * np.cos(front_rotation[1]),
                                 position_sign * legs_separation[0] / 2 * np.sin(central_rotation[0])]
                self.mesh = Cuboid(self.front_legs_size[1], self.front_legs_size[0], self.front_legs_size[2],
                                   position=mesh_position, rotation=mesh_rotation, rotation_order="XZY")
            if number_of_legs[0] == 3:
                if i < 2:
                    mesh_rotation = [front_rotation[0], central_rotation[0], rotation_sign * front_rotation[1]]
                    mesh_position = [
                        position_sign * legs_separation[0] / 2 * np.cos(central_rotation[0]) + legs_separation[
                            2] / 2 * np.sin(central_rotation[0]),
                        -front_legs_size[1] / 2 * np.cos(front_rotation[0]) * np.cos(front_rotation[1]),
                        -position_sign * legs_separation[0] / 2 * np.sin(central_rotation[0]) + legs_separation[
                            2] / 2 * np.cos(central_rotation[0])]
                    self.mesh = Cuboid(self.front_legs_size[1], self.front_legs_size[0], self.front_legs_size[2],
                                       position=mesh_position, rotation=mesh_rotation, rotation_order="XZY")
                else:
                    mesh_rotation = [rear_rotation[0], central_rotation[0], rear_rotation[1]]
                    mesh_position = [position_sign * legs_separation[2] / 2 * np.sin(central_rotation[0]),
                                     -rear_legs_size[1] / 2 * np.cos(rear_rotation[0]) * np.cos(rear_rotation[1]),
                                     position_sign * legs_separation[2] / 2 * np.cos(central_rotation[0])]
                    self.mesh = Cuboid(self.rear_legs_size[1], self.rear_legs_size[0], self.rear_legs_size[2],
                                       position=mesh_position, rotation=mesh_rotation, rotation_order="XZY")
            if number_of_legs[0] == 4:
                if i < 2:
                    mesh_rotation = [front_rotation[0], central_rotation[0], rotation_sign * front_rotation[1]]
                    mesh_position = [
                        position_sign * legs_separation[0] / 2 * np.cos(central_rotation[0]) + legs_separation[
                            2] / 2 * np.sin(central_rotation[0]),
                        -front_legs_size[1] / 2 * np.cos(front_rotation[0]) * np.cos(front_rotation[1]),
                        -position_sign * legs_separation[0] / 2 * np.sin(central_rotation[0]) + legs_separation[
                            2] / 2 * np.cos(central_rotation[0])]
                    self.mesh = Cuboid(self.front_legs_size[1], self.front_legs_size[0], self.front_legs_size[2],
                                       position=mesh_position, rotation=mesh_rotation, rotation_order="XZY")
                else:
                    mesh_rotation = [rear_rotation[0], central_rotation[0], rotation_sign * rear_rotation[1]]
                    mesh_position = [
                        position_sign * legs_separation[1] / 2 * np.cos(central_rotation[0]) - legs_separation[
                            2] / 2 * np.sin(central_rotation[0]),
                        -rear_legs_size[1] / 2 * np.cos(rear_rotation[0]) * np.cos(rear_rotation[1]),
                        -position_sign * legs_separation[1] / 2 * np.sin(central_rotation[0]) - legs_separation[
                            2] / 2 * np.cos(central_rotation[0])]
                    self.mesh = Cuboid(self.rear_legs_size[1], self.rear_legs_size[0], self.rear_legs_size[2],
                                       position=mesh_position, rotation=mesh_rotation, rotation_order="XZY")

            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Leg'


class C_shaped_office_leg(ConceptTemplate):
    def __init__(self, vertical_leg_size, horizontal_z_leg_size, horizontal_x_leg_size,
                 vertical_leg_separation, vertical_leg_rotation, horizontal_leg_rotation,
                 position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        vertical_leg_rotation = [x / 180 * np.pi for x in vertical_leg_rotation]
        horizontal_leg_rotation = [x / 180 * np.pi for x in horizontal_leg_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.vertical_leg_size = vertical_leg_size
        self.horizontal_z_leg_size = horizontal_z_leg_size
        self.horizontal_x_leg_size = horizontal_x_leg_size
        self.vertical_leg_separation = vertical_leg_separation
        self.vertical_leg_rotation = vertical_leg_rotation
        self.horizontal_leg_rotation = horizontal_leg_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(5):
            rotation_sign = -1 if (i % 2 == 1) else 1
            position_sign = -1 if (i % 2 == 1) else 1
            if i < 2:
                mesh_rotation = [vertical_leg_rotation[0], 0, rotation_sign * vertical_leg_rotation[1]]
                mesh_position = [position_sign * vertical_leg_separation[0] / 2,
                                 -vertical_leg_size[1] / 2 * np.cos(vertical_leg_rotation[0]) * np.cos(
                                     vertical_leg_rotation[1]), 0]
                self.mesh = Cuboid(self.vertical_leg_size[1], self.vertical_leg_size[0], self.vertical_leg_size[2],
                                   position=mesh_position, rotation=mesh_rotation, rotation_order="ZXY")
            elif i < 4:
                mesh_rotation = [horizontal_leg_rotation[0], rotation_sign * horizontal_leg_rotation[1], 0]
                mesh_position = [position_sign * (vertical_leg_separation[0] / 2 + vertical_leg_size[1] / 2 * np.sin(
                    vertical_leg_rotation[1]) - horizontal_z_leg_size[1] / 2 * np.sin(horizontal_leg_rotation[1])),
                                 -(vertical_leg_size[1] - horizontal_z_leg_size[0] / 2) * np.cos(
                                     vertical_leg_rotation[0]) * np.cos(
                                     vertical_leg_rotation[1]) + (
                                         horizontal_z_leg_size[1] + horizontal_z_leg_size[0]) / 2 * np.sin(
                                     horizontal_leg_rotation[0]),
                                 -vertical_leg_size[1] / 2 * np.cos(vertical_leg_rotation[1]) * np.sin(
                                     vertical_leg_rotation[0]) - (
                                         vertical_leg_size[2] + horizontal_z_leg_size[1]) / 2 * np.cos(
                                     horizontal_leg_rotation[1]) * np.cos(horizontal_leg_rotation[0])]
                self.mesh = Cuboid(self.horizontal_z_leg_size[0], self.vertical_leg_size[0],
                                   self.horizontal_z_leg_size[1],
                                   position=mesh_position, rotation=mesh_rotation, rotation_order="YXZ")
            else:
                mesh_rotation = [horizontal_leg_rotation[0], 0, 0]
                mesh_position = [0,
                                 -(vertical_leg_size[1] - horizontal_z_leg_size[0] / 2) * np.cos(
                                     vertical_leg_rotation[0]) * np.cos(
                                     vertical_leg_rotation[1]) + (
                                         horizontal_z_leg_size[1] * 2 + horizontal_z_leg_size[0] +
                                         horizontal_x_leg_size[0]) / 2 * np.sin(horizontal_leg_rotation[0]),
                                 -vertical_leg_size[1] / 2 * np.cos(vertical_leg_rotation[1]) * np.sin(
                                     vertical_leg_rotation[0]) - (vertical_leg_size[2] + horizontal_z_leg_size[1] * 2 +
                                                                  horizontal_x_leg_size[0]) / 2 * np.cos(
                                     horizontal_leg_rotation[1]) * np.cos(horizontal_leg_rotation[0])]
                self.mesh = Cuboid(self.horizontal_z_leg_size[0],
                                   vertical_leg_separation[0] + vertical_leg_size[0] + vertical_leg_size[1] * np.sin(
                                       vertical_leg_rotation[1]) - 2 * (
                                           horizontal_z_leg_size[1] + vertical_leg_size[2] / 2 +
                                           horizontal_x_leg_size[0] / 2) * np.sin(horizontal_leg_rotation[1]),
                                   self.horizontal_x_leg_size[0],
                                   position=mesh_position, rotation=mesh_rotation)

            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Leg'


class Star_leg(ConceptTemplate):
    def __init__(self, vertical_sizes, sub_sizes, sub_central_offset,
                 tilt_angle, central_rotation, horizontal_rotation, number_of_sub_legs,
                 position=[0, 0, 0], rotation=[0, 0, 0]):
        
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        tilt_angle = [x / 180 * np.pi for x in tilt_angle]
        central_rotation = [x / 180 * np.pi for x in central_rotation]
        horizontal_rotation = [x / 180 * np.pi for x in horizontal_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.vertical_sizes = vertical_sizes
        self.sub_sizes = sub_sizes
        self.sub_central_offset = sub_central_offset
        self.tilt_angle = tilt_angle
        self.central_rotation = central_rotation
        self.horizontal_rotation = horizontal_rotation
        self.number_of_sub_legs = number_of_sub_legs

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(1 + number_of_sub_legs[0]):
            if i == 0:
                mesh_rotation = [horizontal_rotation[0], central_rotation[0], 0]
                mesh_position = [0, 0, 0]
                self.mesh = Cylinder(vertical_sizes[1], vertical_sizes[0], position=mesh_position,
                                     rotation=mesh_rotation)
            else:
                sub_rotation = np.pi / number_of_sub_legs[0] * (i - 1) * 2
                mesh_position = [
                    -(sub_sizes[2] / 2 * np.cos(tilt_angle[0]) * np.sin(sub_rotation)) * np.cos(central_rotation[0]) + (
                            sub_sizes[2] / 2 * np.cos(tilt_angle[0]) * np.cos(sub_rotation)) * np.sin(
                        central_rotation[0]),
                    (-vertical_sizes[1] / 2 + sub_sizes[1] / 2 + sub_central_offset[0]) * np.cos(horizontal_rotation[0]) - (
                            (sub_sizes[2] / 2 * np.cos(tilt_angle[0]) * np.cos(sub_rotation)) * np.cos(
                        central_rotation[0]) + (
                                    sub_sizes[2] / 2 * np.cos(tilt_angle[0]) * np.sin(sub_rotation)) * np.sin(
                        central_rotation[0])) * np.sin(horizontal_rotation[0]),
                    ((sub_sizes[2] / 2 * np.cos(tilt_angle[0]) * np.cos(sub_rotation)) * np.cos(central_rotation[0]) + (
                            sub_sizes[2] / 2 * np.cos(tilt_angle[0]) * np.sin(sub_rotation)) * np.sin(
                        central_rotation[0])) * np.cos(horizontal_rotation[0]) + (
                            -vertical_sizes[1] + sub_sizes[1] / 2 + sub_central_offset[0]) * np.sin(
                        horizontal_rotation[0]) + vertical_sizes[1] / 2 * np.sin(horizontal_rotation[0])]
                self.mesh = Cuboid(sub_sizes[1], sub_sizes[0], sub_sizes[2])

                tilt_mat = np.array(get_rodrigues_matrix([1, 0, 0], tilt_angle[0]))
                self.mesh.vertices = np.matmul(self.mesh.vertices, tilt_mat.T)
                sub_mat = np.array(get_rodrigues_matrix([0, 1, 0], -sub_rotation))
                self.mesh.vertices = np.matmul(self.mesh.vertices, sub_mat.T)

                cen_mat4 = np.array(get_rodrigues_matrix([0, 1, 0], central_rotation[0]))
                self.mesh.vertices = np.matmul(self.mesh.vertices, cen_mat4.T)

                h_mat4 = np.array(get_rodrigues_matrix([1, 0, 0], horizontal_rotation[0]))
                self.mesh.vertices = np.matmul(self.mesh.vertices, h_mat4.T)

                self.mesh.vertices = apply_transformation(self.mesh.vertices, mesh_position, [0, 0, 0])


            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Leg'


class Regular_leg_with_splat(ConceptTemplate):
    def __init__(self, front_legs_size, rear_legs_size, legs_separation,
                 central_rotation, front_rotation, rear_rotation,
                 front_rear_bridging_bars_sizes, left_right_bridging_bars_sizes, front_rear_bridging_bars_offset,
                 left_right_bridging_bars_offset, bridging_bars_existance, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        central_rotation = [x / 180 * np.pi for x in central_rotation]
        front_rotation = [x / 180 * np.pi for x in front_rotation]
        rear_rotation = [x / 180 * np.pi for x in rear_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.front_legs_size = front_legs_size
        self.rear_legs_size = rear_legs_size
        self.legs_separation = legs_separation
        self.central_rotation = central_rotation
        self.front_rotation = front_rotation
        self.rear_rotation = rear_rotation
        self.front_rear_bridging_bars_sizes = front_rear_bridging_bars_sizes
        self.left_right_bridging_bars_sizes = left_right_bridging_bars_sizes
        self.front_rear_bridging_bars_offset = front_rear_bridging_bars_offset
        self.left_right_bridging_bars_offset = left_right_bridging_bars_offset
        self.bridging_bars_existance = bridging_bars_existance

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(4):
            rotation_sign = -1 if (i % 2 == 1) else 1
            position_sign = -1 if (i % 2 == 0) else 1
            if i < 2:
                mesh_rotation = [front_rotation[0], central_rotation[0], rotation_sign * front_rotation[1]]
                mesh_position = [
                    position_sign * legs_separation[0] / 2 * np.cos(central_rotation[0]) + legs_separation[
                        2] / 2 * np.sin(central_rotation[0]),
                    -front_legs_size[1] / 2 * np.cos(front_rotation[0]) * np.cos(front_rotation[1]),
                    -position_sign * legs_separation[0] / 2 * np.sin(central_rotation[0]) + legs_separation[
                        2] / 2 * np.cos(central_rotation[0])]
                self.mesh = Cuboid(self.front_legs_size[1], self.front_legs_size[0], self.front_legs_size[2],
                                   position=mesh_position, rotation=mesh_rotation, rotation_order="XZY")
            else:
                mesh_rotation = [rear_rotation[0], central_rotation[0], rotation_sign * rear_rotation[1]]
                mesh_position = [
                        position_sign * legs_separation[1] / 2 * np.cos(central_rotation[0]) - legs_separation[
                            2] / 2 * np.sin(central_rotation[0]),
                        -rear_legs_size[1] / 2 * np.cos(rear_rotation[0]) * np.cos(rear_rotation[1]),
                        -position_sign * legs_separation[1] / 2 * np.sin(central_rotation[0]) - legs_separation[
                            2] / 2 * np.cos(central_rotation[0])]
                self.mesh = Cuboid(self.rear_legs_size[1], self.rear_legs_size[0], self.rear_legs_size[2],
                                   position=mesh_position, rotation=mesh_rotation, rotation_order="XZY")

            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        if bridging_bars_existance[0] == 1:
            mesh_rotation = [front_rotation[0], central_rotation[0], 0]
            mesh_position = [
                (legs_separation[2] / 2 + front_rear_bridging_bars_offset[0] * np.sin(front_rotation[0])) * np.sin(
                    central_rotation[0]),
                -(front_legs_size[1] / 2 - front_rear_bridging_bars_offset[0]) * np.cos(front_rotation[0]) * np.cos(
                    front_rotation[1]),
                (legs_separation[2] / 2 + front_rear_bridging_bars_offset[0] * np.sin(front_rotation[0])) * np.cos(
                    central_rotation[0])]
            self.mesh = Cuboid(front_rear_bridging_bars_sizes[0],
                               legs_separation[0] - front_legs_size[0] + front_rear_bridging_bars_offset[0] * np.sin(
                                   front_rotation[1]) * 2,
                               front_rear_bridging_bars_sizes[1],
                               position=mesh_position, rotation=mesh_rotation, rotation_order="XZY")
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        if bridging_bars_existance[1] == 1:
            mesh_rotation = [rear_rotation[0], central_rotation[0], 0]
            mesh_position = [
                (-legs_separation[2] / 2 + front_rear_bridging_bars_offset[1] * np.sin(rear_rotation[0])) * np.sin(
                    central_rotation[0]),
                -(rear_legs_size[1] / 2 - front_rear_bridging_bars_offset[1]) * np.cos(rear_rotation[0]) * np.cos(
                    rear_rotation[1]),
                (-legs_separation[2] / 2 + front_rear_bridging_bars_offset[1] * np.sin(rear_rotation[0])) * np.cos(
                    central_rotation[0])]
            self.mesh = Cuboid(front_rear_bridging_bars_sizes[0],
                               legs_separation[1] - rear_legs_size[0] + front_rear_bridging_bars_offset[1] * np.sin(
                                   front_rotation[1]) * 2,
                               front_rear_bridging_bars_sizes[1],
                               position=mesh_position, rotation=mesh_rotation, rotation_order="XZY")
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        if bridging_bars_existance[2] == 1:

            fr_an_x = - legs_separation[0] / 2 - left_right_bridging_bars_offset[0] * np.cos(
                front_rotation[0]) * np.sin(front_rotation[1])
            fr_an_y = -(front_legs_size[1] / 2 - left_right_bridging_bars_offset[0]) * np.cos(
                front_rotation[0]) * np.cos(front_rotation[1])
            fr_an_z = legs_separation[2] / 2 + left_right_bridging_bars_offset[0] * np.sin(
                front_rotation[0])

            re_an_x = - legs_separation[1] / 2 - left_right_bridging_bars_offset[1] * np.cos(
                rear_rotation[0]) * np.sin(rear_rotation[1])
            re_an_y = -(rear_legs_size[1] / 2 - left_right_bridging_bars_offset[1]) * np.cos(
                rear_rotation[0]) * np.cos(rear_rotation[1])
            re_an_z = - legs_separation[2] / 2 + left_right_bridging_bars_offset[1] * np.sin(
                rear_rotation[0])

            diff_x, diff_y, diff_z = fr_an_x - re_an_x, fr_an_y - re_an_y, fr_an_z - re_an_z
            diff_norm = np.sqrt(diff_x ** 2 + diff_y ** 2 + diff_z ** 2)

            if diff_norm > 0:
                mesh_rotation = [-np.arcsin(diff_y / diff_norm), np.arctan(diff_x / np.sign(diff_z) / (np.abs(diff_z) + 1e-7)) + central_rotation[0], 0]
                mesh_position = [
                    (fr_an_x + re_an_x) / 2 * np.cos(central_rotation[0]) + (fr_an_z + re_an_z) / 2 * np.sin(
                        central_rotation[0]),
                    (fr_an_y + re_an_y) / 2,
                    (fr_an_z + re_an_z) / 2 * np.cos(central_rotation[0]) - (fr_an_x + re_an_x) / 2 * np.sin(
                        central_rotation[0])
                ]

                self.mesh = Cuboid(
                    left_right_bridging_bars_sizes[1],
                    left_right_bridging_bars_sizes[0],
                    diff_norm - (front_legs_size[2] + rear_legs_size[2]) / 2,
                    position=mesh_position,
                    rotation=mesh_rotation, rotation_order="XZY"
                )
                vertices_list.append(self.mesh.vertices)
                faces_list.append(self.mesh.faces + total_num_vertices)
                total_num_vertices += len(self.mesh.vertices)

        if bridging_bars_existance[3] == 1:

            fr_an_x = legs_separation[0] / 2 + left_right_bridging_bars_offset[0] * np.cos(
                front_rotation[0]) * np.sin(front_rotation[1])
            fr_an_y = -(front_legs_size[1] / 2 - left_right_bridging_bars_offset[0]) * np.cos(
                front_rotation[0]) * np.cos(front_rotation[1])
            fr_an_z = legs_separation[2] / 2 + left_right_bridging_bars_offset[0] * np.sin(
                front_rotation[0])

            re_an_x = legs_separation[1] / 2 + left_right_bridging_bars_offset[1] * np.cos(
                rear_rotation[0]) * np.sin(rear_rotation[1])
            re_an_y = -(front_legs_size[1] / 2 - left_right_bridging_bars_offset[1]) * np.cos(
                rear_rotation[0]) * np.cos(rear_rotation[1])
            re_an_z = - legs_separation[2] / 2 + left_right_bridging_bars_offset[1] * np.sin(
                rear_rotation[0])

            diff_x, diff_y, diff_z = fr_an_x - re_an_x, fr_an_y - re_an_y, fr_an_z - re_an_z
            diff_norm = np.sqrt(diff_x ** 2 + diff_y ** 2 + diff_z ** 2)

            if diff_norm > 0:

                mesh_rotation = [-np.arcsin(diff_y / diff_norm), np.arctan(diff_x / np.sign(diff_z) / (np.abs(diff_z) + 1e-7)) + central_rotation[0], 0]
                mesh_position = [
                    (fr_an_x + re_an_x) / 2 * np.cos(central_rotation[0]) + (fr_an_z + re_an_z) / 2 * np.sin(
                        central_rotation[0]),
                    (fr_an_y + re_an_y) / 2,
                    (fr_an_z + re_an_z) / 2 * np.cos(central_rotation[0]) - (fr_an_x + re_an_x) / 2 * np.sin(
                        central_rotation[0])
                ]

                self.mesh = Cuboid(
                    left_right_bridging_bars_sizes[1],
                    left_right_bridging_bars_sizes[0],
                    diff_norm - (front_legs_size[2] + rear_legs_size[2]) / 2,
                    position=mesh_position,
                    rotation=mesh_rotation, rotation_order="XZY"
                )
                vertices_list.append(self.mesh.vertices)
                faces_list.append(self.mesh.faces + total_num_vertices)
                total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Leg'


class Barstool_leg(ConceptTemplate):
    def __init__(self, vertical_sizes, bottom_sizes, horizontal_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):
        
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        horizontal_rotation = [x / 180 * np.pi for x in horizontal_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.vertical_sizes = vertical_sizes
        self.bottom_sizes = bottom_sizes
        self.horizontal_rotation = horizontal_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_rotation = [horizontal_rotation[0], 0, 0]
        mesh_position = [0, -vertical_sizes[1] / 2 * np.cos(horizontal_rotation[0]), 0]
        self.support_mesh = Cylinder(vertical_sizes[1], vertical_sizes[0], position=mesh_position,
                                     rotation=mesh_rotation)
        vertices_list.append(self.support_mesh.vertices)
        faces_list.append(self.support_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.support_mesh.vertices)

        mesh_rotation = [horizontal_rotation[0], 0, 0]
        mesh_position = [0, (-vertical_sizes[1] - bottom_sizes[1] / 2) * np.cos(horizontal_rotation[0]),
                         (vertical_sizes[1] - bottom_sizes[1]) / 2 * np.sin(horizontal_rotation[0])]
        self.bottom_mesh = Cylinder(bottom_sizes[1], bottom_sizes[0], position=mesh_position,
                                    rotation=mesh_rotation)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Leg'


class Solid_armrest(ConceptTemplate):
    def __init__(self, size, armrest_separation, armrest_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        armrest_rotation = [x / 180 * np.pi for x in armrest_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.armrest_separation = armrest_separation
        self.armrest_rotation = armrest_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(2):
            flag = 1 if i == 0 else -1
            mesh_rotation = [armrest_rotation[0], 0, flag * armrest_rotation[1]]
            mesh_position = [-flag * armrest_separation[0] / 2,
                             size[1] / 2 * np.cos(armrest_rotation[1]) * np.cos(armrest_rotation[0]) - size[
                                 2] / 2 * np.sin(armrest_rotation[0]), 0]
            self.mesh = Cuboid(size[1], size[0], size[2], position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Armrest'


class Office_armrest(ConceptTemplate):
    def __init__(self, horizontal_support_sizes, vertical_support_sizes, supports_contact_offset,
                 vertical_support_rotation, horizontal_support_rotation, armrest_separation, armrest_rotation,
                 position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        vertical_support_rotation = [x / 180 * np.pi for x in vertical_support_rotation]
        horizontal_support_rotation = [x / 180 * np.pi for x in horizontal_support_rotation]
        armrest_rotation = [x / 180 * np.pi for x in armrest_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.horizontal_support_sizes = horizontal_support_sizes
        self.vertical_support_sizes = vertical_support_sizes
        self.supports_contact_offset = supports_contact_offset
        self.vertical_support_rotation = vertical_support_rotation
        self.horizontal_support_rotation = horizontal_support_rotation
        self.armrest_separation = armrest_separation
        self.armrest_rotation = armrest_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(4):
            flag = 1 if (i % 2 == 0) else -1
            if i < 2:
                mesh_rotation = [horizontal_support_rotation[0], 0, flag * armrest_rotation[0]]
                mesh_position = [-flag * (armrest_separation[0] / 2 + (
                        vertical_support_sizes[1] * np.cos(vertical_support_rotation[0]) +
                        horizontal_support_sizes[1]) / 2 * np.sin(armrest_rotation[0])),
                                 (vertical_support_sizes[1] * np.cos(vertical_support_rotation[0]) * 2 +
                                  horizontal_support_sizes[1]) / 2 * np.cos(armrest_rotation[0]) + (
                                         supports_contact_offset[0] - vertical_support_sizes[1] / 2 * np.sin(
                                     vertical_support_rotation[0])) * np.tan(horizontal_support_rotation[0]),
                                 0]
                self.mesh = Cuboid(horizontal_support_sizes[1], horizontal_support_sizes[0],
                                   horizontal_support_sizes[2],
                                   position=mesh_position, rotation=mesh_rotation, rotation_order="YXZ")
            else:
                mesh_rotation = [vertical_support_rotation[0], 0, flag * armrest_rotation[0]]
                mesh_position = [-flag * armrest_separation[0] / 2,
                                 vertical_support_sizes[1] / 2 * np.cos(vertical_support_rotation[0]) * np.cos(
                                     armrest_rotation[0]),
                                 supports_contact_offset[0]]
                self.mesh = Cuboid(vertical_support_sizes[1], vertical_support_sizes[0],
                                   vertical_support_sizes[2],
                                   position=mesh_position, rotation=mesh_rotation, rotation_order="YXZ")

            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Armrest'
