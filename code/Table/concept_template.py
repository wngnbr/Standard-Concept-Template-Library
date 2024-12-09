import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation, get_rodrigues_matrix
from knowledge_utils import *
import trimesh


class Regular_desktop(ConceptTemplate):
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

        self.desktop_mesh = Cuboid(self.size[1], self.size[0], self.size[2])
        vertices_list.append(self.desktop_mesh.vertices)
        faces_list.append(self.desktop_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.desktop_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Cylindrical_desktop(ConceptTemplate):
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

        self.desktop_mesh = Cylinder(self.height, self.radius)
        vertices_list.append(self.desktop_mesh.vertices)
        faces_list.append(self.desktop_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.desktop_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class L_type_desktop(ConceptTemplate):
    def __init__(self, horizontal_size, vertical_size, desktop_angle, position=[0, 0, 0], rotation=[0, 0, 0]):
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        desktop_angle = [x / 180 * np.pi for x in desktop_angle]

        # Record Parameters
        self.horizontal_size = horizontal_size
        self.vertical_size = vertical_size
        self.desktop_angle = desktop_angle
        self.position = position
        self.rotation = rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.base_mesh = Cuboid(self.horizontal_size[1], self.horizontal_size[0], self.horizontal_size[2])
        vertices_list.append(self.base_mesh.vertices)
        faces_list.append(self.base_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.base_mesh.vertices)

        other_mesh_rotation = [0, desktop_angle[0], 0]
        other_mesh_position = [-horizontal_size[0] / 2 - np.sqrt(vertical_size[0] ** 2) *
                               np.sin(desktop_angle[0]),
                               0,
                               -horizontal_size[1] / 2 - np.sqrt(vertical_size[0] ** 2) *
                               np.cos(desktop_angle[0])]
        self.other_mesh = Cuboid(self.horizontal_size[1], self.horizontal_size[0], self.horizontal_size[2],
                                 position=other_mesh_position, rotation=other_mesh_rotation)
        vertices_list.append(self.other_mesh.vertices)
        faces_list.append(self.other_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.other_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Regular_leg(ConceptTemplate):
    def __init__(self, front_legs_size, rear_legs_size, number_of_legs, legs_separation,
                 central_rotation, front_rotation, rear_rotation, symmetry_mode, additional_legs_params,
                 position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        central_rotation = [x / 180 * np.pi for x in central_rotation]
        front_rotation = [x / 180 * np.pi for x in front_rotation]
        rear_rotation = [x / 180 * np.pi for x in rear_rotation]
        additional_legs_params = [x / 180 * np.pi if ((i > 1) and ((i - 1) % 9 in [6, 7, 8])) else x for i, x in
                                  enumerate(additional_legs_params)]
        super().__init__(position, rotation)

        # Record Parameters
        self.front_legs_size = front_legs_size
        self.rear_legs_size = rear_legs_size
        self.number_of_legs = number_of_legs
        self.legs_separation = legs_separation
        self.central_rotation = central_rotation
        self.front_rotation = front_rotation
        self.rear_rotation = rear_rotation
        self.symmetry_mode = symmetry_mode
        self.number_of_additional_legs = additional_legs_params[0]
        self.additional_legs_attributes = additional_legs_params[1:]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(number_of_legs[0]):
            if number_of_legs[0] == 1:
                mesh_rotation = [front_rotation[0], central_rotation[0], 0]
                mesh_position = [0, -front_legs_size[1] * np.cos(front_rotation[0]) / 2, 0]
                self.mesh = Cuboid(self.front_legs_size[1], self.front_legs_size[0],
                                   self.front_legs_size[2],
                                   position=mesh_position, rotation=mesh_rotation)

            if number_of_legs[0] == 2:
                rotation_sign = 1 if (i == 0 or self.symmetry_mode == 0) else -1
                position_sign = 1 if i == 1 else -1
                mesh_rotation = [front_rotation[0], central_rotation[0], rotation_sign * front_rotation[1]]
                mesh_position = [position_sign * legs_separation[0] / 2 * np.cos(central_rotation[0]),
                                 -front_legs_size[1] / 2 * np.cos(front_rotation[0]),
                                 position_sign * legs_separation[0] / 2 * np.sin(central_rotation[0])]
                self.mesh = Cuboid(self.front_legs_size[1], self.front_legs_size[0],
                                   self.front_legs_size[2],
                                   position=mesh_position, rotation=mesh_rotation)

            if number_of_legs[0] == 3:
                rotation_sign = 1 if (i == 0 or self.symmetry_mode == 0) else -1
                position_sign = 1 if i == 1 else -1
                if i < 2:
                    mesh_rotation = [front_rotation[0], central_rotation[0], rotation_sign * front_rotation[1]]
                    mesh_position = [
                        position_sign * legs_separation[0] / 2 * np.cos(central_rotation[0]) + legs_separation[
                            2] / 2 * np.sin(central_rotation[0]),
                        -front_legs_size[1] / 2 * np.cos(front_rotation[0]) * np.cos(front_rotation[1]),
                        -position_sign * legs_separation[0] / 2 * np.sin(central_rotation[0]) + legs_separation[
                            2] / 2 * np.cos(central_rotation[0])]
                    size = front_legs_size
                else:
                    mesh_rotation = [rear_rotation[0], central_rotation[0], rear_rotation[1]]
                    mesh_position = [-legs_separation[2] / 2 * np.sin(central_rotation[0]),
                                     -rear_legs_size[1] / 2 * np.cos(rear_rotation[0]) * np.cos(rear_rotation[1]),
                                     -legs_separation[2] / 2 * np.cos(central_rotation[0])]
                    size = rear_legs_size
                self.mesh = Cuboid(size[1], size[0], size[2],
                                   position=mesh_position, rotation=mesh_rotation)

            if number_of_legs[0] == 4:
                rotation_sign = -1 if (i % 2 == 1 and symmetry_mode[0] == 1) else 1
                position_sign = 1 if (i % 2 == 1) else -1
                if i < 2:
                    mesh_rotation = [front_rotation[0], central_rotation[0], rotation_sign * front_rotation[1]]
                    mesh_position = [
                        position_sign * legs_separation[0] / 2 * np.cos(central_rotation[0]) + legs_separation[
                            2] / 2 * np.sin(central_rotation[0]),
                        -front_legs_size[1] / 2 * np.cos(front_rotation[0]) * np.cos(front_rotation[1]),
                        -position_sign * legs_separation[0] / 2 * np.sin(central_rotation[0]) + legs_separation[
                            2] / 2 * np.cos(central_rotation[0])]
                    size = front_legs_size
                else:
                    mesh_rotation = [rear_rotation[0], central_rotation[0], rotation_sign * rear_rotation[1]]
                    mesh_position = [
                        position_sign * legs_separation[1] / 2 * np.cos(central_rotation[0]) - legs_separation[
                            2] / 2 * np.sin(central_rotation[0]),
                        -front_legs_size[1] / 2 * np.cos(rear_rotation[0]) * np.cos(rear_rotation[1]),
                        -position_sign * legs_separation[1] / 2 * np.sin(central_rotation[0]) - legs_separation[
                            2] / 2 * np.cos(central_rotation[0])]
                    size = rear_legs_size
                self.mesh = Cuboid(size[1], size[0], size[2],
                                   position=mesh_position, rotation=mesh_rotation)

            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        for i in range(self.number_of_additional_legs):
            mesh_position = [self.additional_legs_attributes[9 * i + 3] - position[0],
                             self.additional_legs_attributes[9 * i + 4] -
                             self.additional_legs_attributes[9 * i + 1] * np.cos(self.additional_legs_attributes[9 * i + 6]) * np.cos(
                                 self.additional_legs_attributes[9 * i + 8]) / 2,
                             self.additional_legs_attributes[9 * i + 5] - position[2]]
            mesh_rotation = [self.additional_legs_attributes[9 * i + 6],
                             self.additional_legs_attributes[9 * i + 7],
                             self.additional_legs_attributes[9 * i + 8]]
            self.additional_mesh = Cuboid(self.additional_legs_attributes[9 * i + 1],
                                          self.additional_legs_attributes[9 * i],
                                          self.additional_legs_attributes[9 * i + 2],
                                          position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.additional_mesh.vertices)
            faces_list.append(self.additional_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.additional_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Leg'


class Regular_with_splat_leg(ConceptTemplate):
    def __init__(self, front_legs_size, rear_legs_size, legs_separation,
                 central_rotation, front_rotation, rear_rotation,
                 front_rear_bridging_bars_sizes, left_right_bridging_bars_sizes,
                 front_rear_bridging_bars_offset, left_right_bridging_bars_offset, bridging_bars_existance,
                 additional_legs_params, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        central_rotation = [x / 180 * np.pi for x in central_rotation]
        front_rotation = [x / 180 * np.pi for x in front_rotation]
        rear_rotation = [x / 180 * np.pi for x in rear_rotation]
        additional_legs_params = [x / 180 * np.pi if ((i > 1) and ((i - 1) % 9 in [6, 7, 8])) else x for i, x in
                                  enumerate(additional_legs_params)]
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
        self.number_of_additional_legs = additional_legs_params[0]
        self.additional_legs_attributes = additional_legs_params[1:]

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
                                   position=mesh_position, rotation=mesh_rotation)
            else:
                mesh_rotation = [rear_rotation[0], central_rotation[0], rotation_sign * rear_rotation[1]]
                mesh_position = [
                    position_sign * legs_separation[1] / 2 * np.cos(central_rotation[0]) - legs_separation[
                        2] / 2 * np.sin(central_rotation[0]),
                    -front_legs_size[1] / 2 * np.cos(front_rotation[0]) * np.cos(front_rotation[1]),
                    -position_sign * legs_separation[1] / 2 * np.sin(central_rotation[0]) - legs_separation[
                        2] / 2 * np.cos(central_rotation[0])]
                self.mesh = Cuboid(self.rear_legs_size[1], self.rear_legs_size[0], self.rear_legs_size[2],
                                   position=mesh_position, rotation=mesh_rotation)

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
                               position=mesh_position, rotation=mesh_rotation)
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
                               position=mesh_position, rotation=mesh_rotation)
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
                mesh_rotation = [-np.arcsin(diff_y / diff_norm),
                                 np.arctan(diff_x / np.sign(diff_z) / (np.abs(diff_z) + 1e-7)) + central_rotation[0], 0]
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
                    rotation=mesh_rotation
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
                mesh_rotation = [-np.arcsin(diff_y / diff_norm),
                                 np.arctan(diff_x / np.sign(diff_z) / (np.abs(diff_z) + 1e-7)) + central_rotation[0], 0]
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
                    rotation=mesh_rotation
                )
                vertices_list.append(self.mesh.vertices)
                faces_list.append(self.mesh.faces + total_num_vertices)
                total_num_vertices += len(self.mesh.vertices)

        for i in range(self.number_of_additional_legs):
            mesh_position = [self.additional_legs_attributes[9 * i + 3] - position[0],
                             self.additional_legs_attributes[9 * i + 4] -
                             self.additional_legs_attributes[9 * i + 1] * np.cos(self.additional_legs_attributes[9 * i + 6]) * np.cos(
                                 self.additional_legs_attributes[9 * i + 8]) / 2,
                             self.additional_legs_attributes[9 * i + 5] - position[2]]
            mesh_rotation = [self.additional_legs_attributes[9 * i + 6],
                             self.additional_legs_attributes[9 * i + 7],
                             self.additional_legs_attributes[9 * i + 8]]
            self.additional_mesh = Cuboid(self.additional_legs_attributes[9 * i + 1],
                                          self.additional_legs_attributes[9 * i],
                                          self.additional_legs_attributes[9 * i + 2],
                                          position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.additional_mesh.vertices)
            faces_list.append(self.additional_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.additional_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Leg'


class Cable_stayed_leg(ConceptTemplate):
    def __init__(self, front_legs_size, rear_legs_size, number_of_legs,
                 connections_size, connections_offset, legs_separation, central_rotation,
                 other_rotation, additional_legs_params, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        central_rotation = [x / 180 * np.pi for x in central_rotation]
        other_rotation = [x / 180 * np.pi for x in other_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.front_legs_size = front_legs_size
        self.rear_legs_size = rear_legs_size
        self.number_of_legs = number_of_legs
        self.connections_size = connections_size
        self.connections_offset = connections_offset
        self.legs_separation = legs_separation
        self.central_rotation = central_rotation[0]
        self.other_rotation = other_rotation
        self.number_of_additional_legs = additional_legs_params[0]
        self.additional_legs_attributes = additional_legs_params[1:]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        if self.number_of_legs[0] == 1:
            mesh_position = [0, -front_legs_size[1] / 2 * np.cos(other_rotation[0]), 0]
            mesh_rotation = [other_rotation[0], central_rotation[0], other_rotation[1]]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        if self.number_of_legs[0] == 4:
            front_position, rear_position = None, None
            for i in range(number_of_legs[0]):
                rotation_sign = 1 if (i == 0) else -1
                position_sign = 1 if (i % 2 == 1) else -1
                if i < 2:
                    mesh_rotation = [other_rotation[0], central_rotation[0], rotation_sign * other_rotation[1]]
                    mesh_position = [
                        position_sign * legs_separation[0] / 2 * np.cos(central_rotation[0]) + legs_separation[
                            2] / 2 * np.sin(central_rotation[0]),
                        -front_legs_size[1] / 2 * np.cos(other_rotation[0]) * np.cos(other_rotation[1]),
                        -position_sign * legs_separation[0] / 2 * np.sin(central_rotation[0]) + legs_separation[
                            2] / 2 * np.cos(central_rotation[0])]
                    if i == 0:
                        front_position = mesh_position
                    size = front_legs_size
                else:
                    mesh_rotation = [-other_rotation[0], central_rotation[0], rotation_sign * other_rotation[1]]
                    mesh_position = [
                        position_sign * legs_separation[1] / 2 * np.cos(central_rotation[0]) - legs_separation[
                            2] / 2 * np.sin(central_rotation[0]),
                        -front_legs_size[1] / 2 * np.cos(-other_rotation[0]) * np.cos(other_rotation[1]),
                        -position_sign * legs_separation[1] / 2 * np.sin(central_rotation[0]) - legs_separation[
                            2] / 2 * np.cos(central_rotation[0])]
                    if i == 3:
                        rear_position = mesh_position
                    size = rear_legs_size
                self.mesh = Cuboid(size[1], size[0], size[2],
                                   position=mesh_position, rotation=mesh_rotation)
                vertices_list.append(self.mesh.vertices)
                faces_list.append(self.mesh.faces + total_num_vertices)
                total_num_vertices += len(self.mesh.vertices)
            
            for i in range(2):
                rot_mult = 1 if i == 0 else -1
                length = np.sqrt((front_position[0] - rear_position[0]) ** 2 + (front_position[1] - rear_position[1]) ** 2 + (front_position[2] - rear_position[2]) ** 2)
                direction_vector = np.array([front_position[0] - rear_position[0], front_position[2] - rear_position[2]])
                rotation_angle = np.arctan2(direction_vector[1], direction_vector[0])
                mesh_rotation = [0, central_rotation[0] + rot_mult * rotation_angle, 0]
                mesh_position = [0, -connections_size[1] / 2, 0]
                self.mesh = Cuboid(connections_size[0], length, connections_size[1],
                                   position=mesh_position, rotation=mesh_rotation)
                vertices_list.append(self.mesh.vertices)
                faces_list.append(self.mesh.faces + total_num_vertices)
                total_num_vertices += len(self.mesh.vertices)

        for i in range(self.number_of_additional_legs):
            mesh_position = [self.additional_legs_attributes[9 * i + 3] - position[0],
                             self.additional_legs_attributes[9 * i + 4] -
                             self.additional_legs_attributes[9 * i + 1] * np.cos(self.additional_legs_attributes[9 * i + 6]) * np.cos(
                                 self.additional_legs_attributes[9 * i + 8]) / 2,
                             self.additional_legs_attributes[9 * i + 5] - position[2]]
            mesh_rotation = [self.additional_legs_attributes[9 * i + 6],
                             self.additional_legs_attributes[9 * i + 7],
                             self.additional_legs_attributes[9 * i + 8]]
            self.additional_mesh = Cuboid(self.additional_legs_attributes[9 * i + 1],
                                          self.additional_legs_attributes[9 * i],
                                          self.additional_legs_attributes[9 * i + 2],
                                          position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.additional_mesh.vertices)
            faces_list.append(self.additional_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.additional_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Leg'


class Star_leg(ConceptTemplate):
    def __init__(self, vertical_size, sub_size, sub_central_offset,
                 tilt_angle, central_rotation, horizontal_rotation, number_of_sub_legs,
                 additional_legs_params, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        tilt_angle = [x / 180 * np.pi for x in tilt_angle]
        central_rotation = [x / 180 * np.pi for x in central_rotation]
        horizontal_rotation = [x / 180 * np.pi for x in horizontal_rotation]
        additional_legs_params = [x / 180 * np.pi if ((i > 1) and ((i - 1) % 9 in [6, 7, 8])) else x for i, x in
                                  enumerate(additional_legs_params)]
        super().__init__(position, rotation)

        # Record Parameters
        self.vertical_size = vertical_size
        self.sub_size = sub_size
        self.sub_central_offset = sub_central_offset
        self.tilt_angle = tilt_angle
        self.central_rotation = central_rotation
        self.horizontal_rotation = horizontal_rotation
        self.number_of_sub_legs = number_of_sub_legs
        self.number_of_additional_legs = additional_legs_params[0]
        self.additional_legs_attributes = additional_legs_params[1:]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(1 + number_of_sub_legs[0]):
            if i == 0:
                mesh_rotation = [horizontal_rotation[0], central_rotation[0], 0]
                mesh_position = [0, -vertical_size[1] / 2 * np.cos(horizontal_rotation[0]), 0]
                self.mesh = Cylinder(vertical_size[1], vertical_size[0], position=mesh_position,
                                     rotation=mesh_rotation)
            else:
                sub_rotation = np.pi / number_of_sub_legs[0] * (i - 1) * 2
                mesh_rotation = [tilt_angle[0] + horizontal_rotation[0], central_rotation[0] - sub_rotation, 0]
                mesh_position = [
                    -(sub_size[2] / 2 * np.cos(tilt_angle[0]) * np.sin(sub_rotation)) * np.cos(central_rotation[0]) + (
                            sub_size[2] / 2 * np.cos(tilt_angle[0]) * np.cos(sub_rotation)) * np.sin(
                        central_rotation[0]),
                    (-vertical_size[1] + sub_size[1] / 2 + sub_central_offset[0]) * np.cos(horizontal_rotation[0]) - (
                            (sub_size[2] / 2 * np.cos(tilt_angle[0]) * np.cos(sub_rotation)) * np.cos(
                        central_rotation[0]) + (
                                    sub_size[2] / 2 * np.cos(tilt_angle[0]) * np.sin(sub_rotation)) * np.sin(
                        central_rotation[0])) * np.sin(horizontal_rotation[0]),
                    ((sub_size[2] / 2 * np.cos(tilt_angle[0]) * np.cos(sub_rotation)) * np.cos(central_rotation[0]) + (
                            sub_size[2] / 2 * np.cos(tilt_angle[0]) * np.sin(sub_rotation)) * np.sin(
                        central_rotation[0])) * np.cos(horizontal_rotation[0]) + (
                            -vertical_size[1] + sub_size[1] / 2 + sub_central_offset[0]) * np.sin(
                        horizontal_rotation[0]) + vertical_size[1] / 2 * np.sin(horizontal_rotation[0])]
                self.mesh = Cuboid(sub_size[1], sub_size[0], sub_size[2], position=mesh_position,
                                   rotation=mesh_rotation)

            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        for i in range(self.number_of_additional_legs):
            mesh_position = [self.additional_legs_attributes[9 * i + 3] - position[0],
                             self.additional_legs_attributes[9 * i + 4] -
                             self.additional_legs_attributes[9 * i + 1] * np.cos(self.additional_legs_attributes[9 * i + 6]) * np.cos(
                                 self.additional_legs_attributes[9 * i + 8]) / 2,
                             self.additional_legs_attributes[9 * i + 5] - position[2]]
            mesh_rotation = [self.additional_legs_attributes[9 * i + 6],
                             self.additional_legs_attributes[9 * i + 7],
                             self.additional_legs_attributes[9 * i + 8]]
            self.additional_mesh = Cuboid(self.additional_legs_attributes[9 * i + 1],
                                          self.additional_legs_attributes[9 * i],
                                          self.additional_legs_attributes[9 * i + 2],
                                          position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.additional_mesh.vertices)
            faces_list.append(self.additional_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.additional_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Leg'


class Bar_cylindrical_leg(ConceptTemplate):
    def __init__(self, vertical_size, bottom_size, horizontal_rotation,
                 additional_legs_params, position=[0, 0, 0], rotation=[0, 0, 0]):
        
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        horizontal_rotation = [x / 180 * np.pi for x in horizontal_rotation]
        additional_legs_params = [x / 180 * np.pi if ((i > 1) and ((i - 1) % 9 in [6, 7, 8])) else x for i, x in
                                  enumerate(additional_legs_params)]
        super().__init__(position, rotation)

        # Record Parameters
        self.vertical_size = vertical_size
        self.bottom_size = bottom_size
        self.horizontal_rotation = horizontal_rotation
        self.number_of_additional_legs = additional_legs_params[0]
        self.additional_legs_attributes = additional_legs_params[1:]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_rotation = [horizontal_rotation[0], 0, 0]
        mesh_position = [0, -vertical_size[1] / 2 * np.cos(horizontal_rotation[0]), 0]
        self.support_mesh = Cylinder(vertical_size[1], vertical_size[0], position=mesh_position,
                                     rotation=mesh_rotation)
        vertices_list.append(self.support_mesh.vertices)
        faces_list.append(self.support_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.support_mesh.vertices)

        mesh_rotation = [horizontal_rotation[0], 0, 0]
        mesh_position = [0, (-vertical_size[1] - bottom_size[1] / 2) * np.cos(horizontal_rotation[0]),
                         (vertical_size[1] - bottom_size[1]) / 2 * np.sin(horizontal_rotation[0])]
        self.bottom_mesh = Cylinder(bottom_size[1], bottom_size[0], position=mesh_position,
                                    rotation=mesh_rotation)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        for i in range(self.number_of_additional_legs):
            mesh_position = [self.additional_legs_attributes[9 * i + 3] - position[0],
                             self.additional_legs_attributes[9 * i + 4] -
                             self.additional_legs_attributes[9 * i + 1] * np.cos(self.additional_legs_attributes[9 * i + 6]) * np.cos(
                                 self.additional_legs_attributes[9 * i + 8]) / 2,
                             self.additional_legs_attributes[9 * i + 5] - position[2]]
            mesh_rotation = [self.additional_legs_attributes[9 * i + 6],
                             self.additional_legs_attributes[9 * i + 7],
                             self.additional_legs_attributes[9 * i + 8]]
            self.additional_mesh = Cuboid(self.additional_legs_attributes[9 * i + 1],
                                          self.additional_legs_attributes[9 * i],
                                          self.additional_legs_attributes[9 * i + 2],
                                          position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.additional_mesh.vertices)
            faces_list.append(self.additional_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.additional_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)
        
        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Leg'


class Bar_cuboid_leg(ConceptTemplate):
    def __init__(self, vertical_size, bottom_size, bottom_rotation, horizontal_rotation,
                 additional_legs_params, position=[0, 0, 0], rotation=[0, 0, 0]):
        
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        bottom_rotation = [x / 180 * np.pi for x in bottom_rotation]
        horizontal_rotation = [x / 180 * np.pi for x in horizontal_rotation]
        additional_legs_params = [x / 180 * np.pi if ((i > 1) and ((i - 1) % 9 in [6, 7, 8])) else x for i, x in
                                  enumerate(additional_legs_params)]
        super().__init__(position, rotation)

        # Record Parameters
        self.vertical_size = vertical_size
        self.bottom_size = bottom_size
        self.bottom_rotation = bottom_rotation
        self.horizontal_rotation = horizontal_rotation
        self.number_of_additional_legs = additional_legs_params[0]
        self.additional_legs_attributes = additional_legs_params[1:]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_rotation = [horizontal_rotation[0], 0, 0]
        mesh_position = [0, -vertical_size[1] / 2 * np.cos(horizontal_rotation[0]), 0]
        self.support_mesh = Cuboid(vertical_size[1], vertical_size[0], vertical_size[2],
                                   position=mesh_position,
                                   rotation=mesh_rotation)
        vertices_list.append(self.support_mesh.vertices)
        faces_list.append(self.support_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.support_mesh.vertices)

        mesh_rotation = [horizontal_rotation[0], bottom_rotation[0], 0]
        mesh_position = [0, (-vertical_size[1] - bottom_size[1] / 2) * np.cos(horizontal_rotation[0]),
                         (vertical_size[1] - bottom_size[1]) / 2 * np.sin(horizontal_rotation[0])]
        self.bottom_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                                  position=mesh_position,
                                  rotation=mesh_rotation)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        for i in range(self.number_of_additional_legs):
            mesh_position = [self.additional_legs_attributes[9 * i + 3] - position[0],
                             self.additional_legs_attributes[9 * i + 4] -
                             self.additional_legs_attributes[9 * i + 1] * np.cos(self.additional_legs_attributes[9 * i + 6]) * np.cos(
                                 self.additional_legs_attributes[9 * i + 8]) / 2,
                             self.additional_legs_attributes[9 * i + 5] - position[2]]
            mesh_rotation = [self.additional_legs_attributes[9 * i + 6],
                             self.additional_legs_attributes[9 * i + 7],
                             self.additional_legs_attributes[9 * i + 8]]
            self.additional_mesh = Cuboid(self.additional_legs_attributes[9 * i + 1],
                                          self.additional_legs_attributes[9 * i],
                                          self.additional_legs_attributes[9 * i + 2],
                                          position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.additional_mesh.vertices)
            faces_list.append(self.additional_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.additional_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Leg'


class Desk_type_leg(ConceptTemplate):
    def __init__(self, vertical_size, horizontal_size, vertical_separation,
                 vertical_rotation, horizontal_rotation, connections_size,
                 number_of_connections, connections_offset, interval_between_connections,
                 additional_legs_params, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        vertical_rotation = [x / 180 * np.pi for x in vertical_rotation]
        horizontal_rotation = [x / 180 * np.pi for x in horizontal_rotation]
        additional_legs_params = [x / 180 * np.pi if ((i > 1) and ((i - 1) % 9 in [6, 7, 8])) else x for i, x in
                                  enumerate(additional_legs_params)]
        super().__init__(position, rotation)

        # Record Parameters
        self.vertical_size = vertical_size
        self.horizontal_size = horizontal_size
        self.vertical_separation = vertical_separation
        self.vertical_rotation = vertical_rotation
        self.horizontal_rotation = horizontal_rotation
        self.connections_size = connections_size
        self.number_of_connections = number_of_connections
        self.connections_offset = connections_offset
        self.interval_between_connections = interval_between_connections
        self.number_of_additional_legs = additional_legs_params[0]
        self.additional_legs_attributes = additional_legs_params[1:]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(4 + number_of_connections[0]):
            rotation_sign = 1 if i == 0 else -1
            position_sign = -1 if (i % 2 == 0) else 1
            pose = apply_transformation([0, -vertical_size[1] / 2, 0], [0, 0, 0],
                                        [vertical_rotation[0], 0, vertical_rotation[1]])
            if i < 2:
                mesh_rotation = [vertical_rotation[0], 0, rotation_sign * vertical_rotation[1]]
                mesh_position = [position_sign * vertical_separation[0] / 2, pose[1], 0]
                self.mesh = Cuboid(vertical_size[1], vertical_size[0], vertical_size[2],
                                   position=mesh_position, rotation=mesh_rotation)
            elif i < 4:
                mesh_rotation = [horizontal_rotation[0], 0, 0]
                mesh_position = [position_sign * (vertical_separation[0] / 2 - pose[0]), 2 * pose[1],
                                 -vertical_size[1] / 2 * np.sin(vertical_rotation[0])]
                self.mesh = Cuboid(horizontal_size[1], horizontal_size[0], horizontal_size[2], position=mesh_position,
                                   rotation=mesh_rotation)
            else:

                mesh_rotation = [vertical_rotation[0], 0, 0]
                mesh_position = [0,
                                 -(vertical_size[1] / 2 - (
                                         (i - 4) * interval_between_connections[0] + connections_offset[
                                     0]) * np.cos(vertical_rotation[0]) * np.cos(vertical_rotation[1])),
                                 ((i - 4) * interval_between_connections[0] + connections_offset[0]) * np.sin(
                                     vertical_rotation[0])]
                self.mesh = Cuboid(connections_size[0],
                                   vertical_separation[0] - vertical_size[0] + 2 * (
                                           connections_offset[0] + (i - 4) * interval_between_connections[0]) * np.sin(
                                       vertical_rotation[1]),
                                   connections_size[1],
                                   position=mesh_position, rotation=mesh_rotation)

            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        for i in range(self.number_of_additional_legs):
            mesh_position = [self.additional_legs_attributes[9 * i + 3] - position[0],
                             self.additional_legs_attributes[9 * i + 4] -
                             self.additional_legs_attributes[9 * i + 1] * np.cos(self.additional_legs_attributes[9 * i + 6]) * np.cos(
                                 self.additional_legs_attributes[9 * i + 8]) / 2,
                             self.additional_legs_attributes[9 * i + 5] - position[2]]
            mesh_rotation = [self.additional_legs_attributes[9 * i + 6],
                             self.additional_legs_attributes[9 * i + 7],
                             self.additional_legs_attributes[9 * i + 8]]
            self.additional_mesh = Cuboid(self.additional_legs_attributes[9 * i + 1],
                                          self.additional_legs_attributes[9 * i],
                                          self.additional_legs_attributes[9 * i + 2],
                                          position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.additional_mesh.vertices)
            faces_list.append(self.additional_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.additional_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Leg'


class Regular_sublayer(ConceptTemplate):
    def __init__(self, subs_size, number_of_subs, subs_offset, interval_between_subs,
                 additional_sublayers_params, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        additional_sublayers_params = [x / 180 * np.pi if ((i > 1) and ((i - 1) % 9 in [6, 7, 8])) else x for i, x in
                                       enumerate(additional_sublayers_params)]
        super().__init__(position, rotation)

        # Record Parameters
        self.subs_size = subs_size
        self.number_of_subs = number_of_subs
        self.subs_offset = subs_offset
        self.interval_between_subs = interval_between_subs
        self.number_of_additional_sublayers = additional_sublayers_params[0]
        self.additional_sublayers_attributes = additional_sublayers_params[1:]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(number_of_subs[0]):
            mesh_rotation = [0, 0, 0]
            mesh_position = [0, subs_offset[0] + i * interval_between_subs[0], 0]
            self.mesh = Cuboid(subs_size[1], subs_size[0], subs_size[2],
                               position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        for i in range(self.number_of_additional_sublayers):
            mesh_position = [self.additional_sublayers_attributes[9 * i + 3] - position[0],
                             self.additional_sublayers_attributes[9 * i + 4] - position[1],
                             self.additional_sublayers_attributes[9 * i + 5] - position[2]]
            mesh_rotation = [self.additional_sublayers_attributes[9 * i + 6],
                             self.additional_sublayers_attributes[9 * i + 7],
                             self.additional_sublayers_attributes[9 * i + 8]]
            self.additional_mesh = Cuboid(self.additional_sublayers_attributes[9 * i + 1],
                                          self.additional_sublayers_attributes[9 * i],
                                          self.additional_sublayers_attributes[9 * i + 2],
                                          position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.additional_mesh.vertices)
            faces_list.append(self.additional_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.additional_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Layer'


class Cylindrical_sublayer(ConceptTemplate):
    def __init__(self, subs_size, number_of_subs, subs_offset, interval_between_subs,
                 additional_sublayers_params, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        additional_sublayers_params = [x / 180 * np.pi if ((i > 1) and ((i - 1) % 9 in [6, 7, 8])) else x for i, x in
                                       enumerate(additional_sublayers_params)]
        super().__init__(position, rotation)

        # Record Parameters
        self.subs_size = subs_size
        self.number_of_subs = number_of_subs
        self.subs_offset = subs_offset
        self.interval_between_subs = interval_between_subs
        self.number_of_additional_sublayers = additional_sublayers_params[0]
        self.additional_sublayers_attributes = additional_sublayers_params[1:]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(number_of_subs[0]):
            mesh_rotation = [0, 0, 0]
            mesh_position = [0, subs_offset[0] + i * interval_between_subs[0], 0]
            self.mesh = Cylinder(subs_size[1], subs_size[0],
                                 position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        for i in range(self.number_of_additional_sublayers):
            mesh_position = [self.additional_sublayers_attributes[9 * i + 3],
                             self.additional_sublayers_attributes[9 * i + 4],
                             self.additional_sublayers_attributes[9 * i + 5]]
            mesh_rotation = [self.additional_sublayers_attributes[9 * i + 6],
                             self.additional_sublayers_attributes[9 * i + 7],
                             self.additional_sublayers_attributes[9 * i + 8]]
            self.additional_mesh = Cuboid(self.additional_sublayers_attributes[9 * i + 1],
                                          self.additional_sublayers_attributes[9 * i],
                                          self.additional_sublayers_attributes[9 * i + 2],
                                          position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.additional_mesh.vertices)
            faces_list.append(self.additional_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.additional_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Layer'


class Regular_backboard(ConceptTemplate):
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

        mesh_rotation = [0, 0, 0]
        mesh_position = [0, -size[1] / 2, 0]
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

        self.semantic = 'Board'


class Regular_drawer(ConceptTemplate):
    def __init__(self, number_of_drawer, drawers_params, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        drawers_params = [x / 180 * np.pi if i % 21 in [14] else x for i, x in enumerate(drawers_params)]
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_drawer = number_of_drawer
        self.drawer_size = [drawers_params[i * 21: i * 21 + 3] for i in range(number_of_drawer[0])]
        self.bottom_size = [drawers_params[i * 21 + 3] for i in range(number_of_drawer[0])]
        self.front_size = [drawers_params[i * 21 + 4: i * 21 + 7] for i in range(number_of_drawer[0])]
        self.front_offset = [drawers_params[i * 21 + 7] for i in range(number_of_drawer[0])]
        self.left_right_inner_size = [drawers_params[i * 21 + 8] for i in range(number_of_drawer[0])]
        self.rear_front_inner_size = [drawers_params[i * 21 + 9] for i in range(number_of_drawer[0])]
        self.number_of_handle = [drawers_params[i * 21 + 10] for i in range(number_of_drawer[0])]
        self.handle_sizes = [drawers_params[i * 21 + 11: i * 21 + 14] for i in range(number_of_drawer[0])]
        self.handle_rotation = [drawers_params[i * 21 + 14] for i in range(number_of_drawer[0])]
        self.handle_offset = [drawers_params[i * 21 + 15: i * 21 + 17] for i in range(number_of_drawer[0])]
        self.handle_separation = [drawers_params[i * 21 + 17] for i in range(number_of_drawer[0])]
        self.drawer_offset = [drawers_params[i * 21 + 18: i * 21 + 21] for i in range(number_of_drawer[0])]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for drawer_idx in range(number_of_drawer[0]):
            for mesh_idx in range(6 + self.number_of_handle[drawer_idx]):
                if mesh_idx < 2:
                    position_sign = -1 if mesh_idx == 0 else 1
                    mesh_position = [position_sign * (
                            self.drawer_size[drawer_idx][0] - self.left_right_inner_size[drawer_idx]) / 2 +
                                     self.drawer_offset[drawer_idx][0],
                                     -self.drawer_size[drawer_idx][1] / 2 + self.drawer_offset[drawer_idx][1],
                                     self.drawer_offset[drawer_idx][2]]
                    self.mesh = Cuboid(self.drawer_size[drawer_idx][1],
                                       self.left_right_inner_size[drawer_idx],
                                       self.drawer_size[drawer_idx][2],
                                       position=mesh_position)
                elif mesh_idx < 4:
                    position_sign = -1 if mesh_idx == 3 else 1
                    mesh_position = [self.drawer_offset[drawer_idx][0],
                                     -self.drawer_size[drawer_idx][1] / 2 + self.drawer_offset[drawer_idx][1],
                                     position_sign * (
                                             self.drawer_size[drawer_idx][2] - self.rear_front_inner_size[drawer_idx]) / 2 + self.drawer_offset[drawer_idx][2]]
                    self.mesh = Cuboid(self.drawer_size[drawer_idx][1],
                                       self.drawer_size[drawer_idx][0] - 2 * self.left_right_inner_size[drawer_idx],
                                       self.rear_front_inner_size[drawer_idx],
                                       position=mesh_position)
                elif mesh_idx == 4:
                    mesh_position = [self.drawer_offset[drawer_idx][0],
                                     -self.drawer_size[drawer_idx][1] + self.drawer_offset[drawer_idx][1] - self.bottom_size[drawer_idx] / 2,
                                     self.drawer_offset[drawer_idx][2]]
                    self.mesh = Cuboid(self.bottom_size[drawer_idx],
                                       self.drawer_size[drawer_idx][0],
                                       self.drawer_size[drawer_idx][2],
                                       position=mesh_position)
                elif mesh_idx == 5:
                    mesh_position = [self.drawer_offset[drawer_idx][0],
                                     -self.drawer_size[drawer_idx][1] / 2 + self.drawer_offset[drawer_idx][1] +
                                     self.front_offset[drawer_idx],
                                     self.drawer_offset[drawer_idx][2] + self.drawer_size[drawer_idx][2] / 2 +
                                     self.front_size[drawer_idx][2] / 2]
                    self.mesh = Cuboid(self.front_size[drawer_idx][1],
                                       self.front_size[drawer_idx][0],
                                       self.front_size[drawer_idx][2],
                                       position=mesh_position)
                else:
                    if self.number_of_handle[drawer_idx] == 2:
                        position_sign = 1 if mesh_idx == 6 else -1
                    else:
                        position_sign = 0
                    mesh_rotation = [0, 0, self.handle_rotation[drawer_idx]]
                    mesh_position = [self.drawer_offset[drawer_idx][0] + self.handle_offset[drawer_idx][0] +
                                     position_sign * self.handle_separation[drawer_idx] / 2,
                                     -self.drawer_size[drawer_idx][1] / 2 + self.drawer_offset[drawer_idx][1] +
                                     self.handle_offset[drawer_idx][0],
                                     self.drawer_offset[drawer_idx][2] + self.drawer_size[drawer_idx][2] / 2 +
                                     self.front_size[drawer_idx][2] + self.front_size[drawer_idx][2] / 2]
                    self.mesh = Cuboid(self.handle_sizes[drawer_idx][1], self.handle_sizes[drawer_idx][0], self.handle_sizes[drawer_idx][2],
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

        self.semantic = 'Drawer'


class Regular_door(ConceptTemplate):
    def __init__(self, number_of_door, doors_params, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        doors_params = [x / 180 * np.pi if i % 13 in [6, 9] else x for i, x in enumerate(doors_params)]
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_door = number_of_door
        self.door_size = [doors_params[i * 13: i * 13 + 3] for i in range(number_of_door[0])]
        self.handle_size = [doors_params[i * 13 + 3: i * 13 + 6] for i in range(number_of_door[0])]
        self.handle_rotation = [doors_params[i * 13 + 6] for i in range(number_of_door[0])]
        self.handle_offset = [doors_params[i * 13 + 7: i * 13 + 9] for i in range(number_of_door[0])]
        self.door_rotation = [doors_params[i * 13 + 9] for i in range(number_of_door[0])]
        self.door_offset = [doors_params[i * 13 + 10: i * 13 + 13] for i in range(number_of_door[0])]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for door_idx in range(self.number_of_door[0]):
            for mesh_idx in range(2):
                if mesh_idx == 0:
                    mesh_rotation = [0, self.door_rotation[door_idx], 0]
                    mesh_position = [self.door_offset[door_idx][0],
                                     self.door_offset[door_idx][1] - self.door_size[door_idx][1] / 2,
                                     self.door_offset[door_idx][2]]
                    self.mesh = Cuboid(self.door_size[door_idx][1], self.door_size[door_idx][0], self.door_size[door_idx][2],
                                       position=mesh_position, rotation=mesh_rotation)
                else:
                    mesh_rotation = [0, self.door_rotation[door_idx], self.handle_rotation[door_idx]]
                    mesh_position = [
                        self.door_offset[door_idx][0] + self.handle_offset[door_idx][0] * np.cos(self.door_rotation[door_idx]) + self.handle_size[door_idx][2] / 2 * np.sin(
                            self.door_rotation[door_idx]),
                        self.door_offset[door_idx][1] - self.door_size[door_idx][1] / 2 + self.handle_offset[door_idx][1],
                        self.door_offset[door_idx][2] - self.handle_offset[door_idx][0] * np.sin(self.door_rotation[door_idx]) + self.handle_size[door_idx][2] / 2 * np.cos(
                            self.door_rotation[door_idx])]
                    self.mesh = Cuboid(self.handle_size[door_idx][1], self.handle_size[door_idx][0], self.handle_size[door_idx][2],
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

        self.semantic = 'Door'


class Regular_cabinet(ConceptTemplate):
    def __init__(self, number_of_top_cabinet, number_of_beneath_cabinet, cab_backs_size,
                 cab_left_right_inner_sizes, cab_up_down_inner_sizes, drawer_inner_sizes,
                 drawer_bottom_size, door_sizes, number_of_layers, layers_sizes, layers_offset,
                 interval_between_layers, cabinets_params, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        cabinets_params = [x / 180 * np.pi if i % 56 in [22, 23, 24, 25, 41, 42, 43, 44] else x for i, x in enumerate(cabinets_params)]
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_top_cabinet = number_of_top_cabinet[0]
        self.number_of_beneath_cabinet = number_of_beneath_cabinet[0]
        self.cab_backs_size = cab_backs_size
        self.cab_left_right_inner_sizes = cab_left_right_inner_sizes
        self.cab_up_down_inner_sizes = cab_up_down_inner_sizes
        self.drawer_inner_sizes = drawer_inner_sizes
        self.drawer_bottom_size = drawer_bottom_size
        self.door_sizes = door_sizes
        self.number_of_layers = number_of_layers
        self.layers_sizes = layers_sizes
        self.layers_offset = layers_offset
        self.interval_between_layers = interval_between_layers
        self.cabinet_size = [cabinets_params[i * 56: i * 56 + 3] for i in range(self.number_of_top_cabinet)
                             ] + [cabinets_params[(i + 2) * 56: (i + 2) * 56 + 3] for i in range(self.number_of_beneath_cabinet)]
        self.type_of_spaces = [
                                  cabinets_params[i * 56 + 3: i * 56 + 7]
                                  for i in range(self.number_of_top_cabinet)
                              ] + [
                                  cabinets_params[(i + 2) * 56 + 3: (i + 2) * 56 + 7]
                                  for i in range(self.number_of_beneath_cabinet)
                              ]

        self.drawer_interval = [
                                   cabinets_params[i * 56 + 7: i * 56 + 11]
                                   for i in range(self.number_of_top_cabinet)
                               ] + [
                                   cabinets_params[(i + 2) * 56 + 6: (i + 2) * 56 + 11]
                                   for i in range(self.number_of_beneath_cabinet)
                               ]

        self.drawer_offset = [
                                 cabinets_params[i * 56 + 11: i * 56 + 15]
                                 for i in range(self.number_of_top_cabinet)
                             ] + [
                                 cabinets_params[(i + 2) * 56 + 11: (i + 2) * 56 + 15]
                                 for i in range(self.number_of_beneath_cabinet)
                             ]

        self.drawer_number_of_handles = [
                                            cabinets_params[i * 56 + 15: i * 56 + 19]
                                            for i in range(self.number_of_top_cabinet)
                                        ] + [
                                            cabinets_params[(i + 2) * 56 + 15: (i + 2) * 56 + 19]
                                            for i in range(self.number_of_beneath_cabinet)
                                        ]

        self.drawer_handles_size = [
                                       cabinets_params[i * 56 + 19: i * 56 + 22]
                                       for i in range(self.number_of_top_cabinet)
                                   ] + [
                                       cabinets_params[(i + 2) * 56 + 19: (i + 2) * 56 + 22]
                                       for i in range(self.number_of_beneath_cabinet)
                                   ]

        self.drawer_handles_rotation = [
                                           cabinets_params[i * 56 + 22: i * 56 + 26]
                                           for i in range(self.number_of_top_cabinet)
                                       ] + [
                                           cabinets_params[(i + 2) * 56 + 22: (i + 2) * 56 + 26]
                                           for i in range(self.number_of_beneath_cabinet)
                                       ]

        self.drawer_handles_separation = [
                                             cabinets_params[i * 56 + 26: i * 56 + 30]
                                             for i in range(self.number_of_top_cabinet)
                                         ] + [
                                             cabinets_params[(i + 2) * 56 + 26: (i + 2) * 56 + 30]
                                             for i in range(self.number_of_beneath_cabinet)
                                         ]

        self.drawer_handles_offsets = [
                                          [
                                              [x_offset, y_offset]
                                              for x_offset, y_offset in zip(
                                              cabinets_params[i * 56 + 30: i * 56 + 34],
                                              cabinets_params[i * 56 + 34: i * 56 + 38]
                                          )
                                          ]
                                          for i in range(self.number_of_top_cabinet)
                                      ] + [
                                          [
                                              [x_offset, y_offset]
                                              for x_offset, y_offset in zip(
                                              cabinets_params[(i + 2) * 56 + 30: (i + 2) * 56 + 34],
                                              cabinets_params[(i + 2) * 56 + 34: (i + 2) * 56 + 38]
                                          )
                                          ]
                                          for i in range(self.number_of_beneath_cabinet)
                                      ]

        self.door_handles_size = [
                                     cabinets_params[i * 56 + 38: i * 56 + 41]
                                     for i in range(self.number_of_top_cabinet)
                                 ] + [
                                     cabinets_params[(i + 2) * 56 + 38: (i + 2) * 56 + 41]
                                     for i in range(self.number_of_beneath_cabinet)
                                 ]

        self.door_handles_rotation = [
                                         cabinets_params[i * 56 + 41: i * 56 + 45]
                                         for i in range(self.number_of_top_cabinet)
                                     ] + [
                                         cabinets_params[(i + 2) * 56 + 41: (i + 2) * 56 + 45]
                                         for i in range(self.number_of_beneath_cabinet)
                                     ]

        self.door_handles_offsets = [
                                        [
                                            [x_offset, y_offset]
                                            for x_offset, y_offset in zip(
                                            cabinets_params[i * 56 + 45: i * 56 + 49],
                                            cabinets_params[i * 56 + 49: i * 56 + 53]
                                        )
                                        ]
                                        for i in range(self.number_of_top_cabinet)
                                    ] + [
                                        [
                                            [x_offset, y_offset]
                                            for x_offset, y_offset in zip(
                                            cabinets_params[(i + 2) * 56 + 45: (i + 2) * 56 + 49],
                                            cabinets_params[(i + 2) * 56 + 49: (i + 2) * 56 + 53]
                                        )
                                        ]
                                        for i in range(self.number_of_beneath_cabinet)
                                    ]

        self.cabinet_offset = [
                                  cabinets_params[i * 56 + 53: i * 56 + 56]
                                  for i in range(self.number_of_top_cabinet)
                              ] + [
                                  cabinets_params[(i + 2) * 56 + 53: (i + 2) * 56 + 56]
                                  for i in range(self.number_of_beneath_cabinet)
                              ]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # top cabinet
        for cabinet_idx in range(self.number_of_top_cabinet):
            cabinet_body_mesh = 5
            layers_mesh = number_of_layers[cabinet_idx]
            cabinet_doors_mesh = sum(2 for type_of_space in self.type_of_spaces[cabinet_idx] if type_of_space == 2)
            cabinet_drawers_mesh = sum(
                (5 + self.drawer_number_of_handles[cabinet_idx][space_idx]) for space_idx, type_of_space in enumerate(self.type_of_spaces[cabinet_idx]) if type_of_space == 1)
            total_mesh_list = [cabinet_body_mesh, layers_mesh, cabinet_doors_mesh, cabinet_drawers_mesh]

            # build cabinet body and layers
            for mesh_idx in range(sum(total_mesh_list[:2])):
                if mesh_idx < sum(total_mesh_list[:1]):
                    position_sign = -1 if mesh_idx % 2 == 0 else 1
                    if mesh_idx < 2:
                        mesh_position = [
                            position_sign * (self.cabinet_size[cabinet_idx][0] - self.cab_left_right_inner_sizes[cabinet_idx]) / 2 + self.cabinet_offset[cabinet_idx][0],
                            self.cabinet_offset[cabinet_idx][1] - self.cabinet_size[cabinet_idx][1] / 2,
                            self.cabinet_offset[cabinet_idx][2]]
                        self.mesh = Cuboid(self.cabinet_size[cabinet_idx][1], self.cab_left_right_inner_sizes[cabinet_idx], self.cabinet_size[cabinet_idx][2],
                                           position=mesh_position)
                    elif mesh_idx < 4:
                        mesh_position = [
                            self.cabinet_offset[cabinet_idx][0],
                            position_sign * (self.cabinet_size[cabinet_idx][1] - self.cab_up_down_inner_sizes[cabinet_idx]) / 2 + self.cabinet_offset[cabinet_idx][1] - self.cabinet_size[cabinet_idx][1] / 2,
                            self.cabinet_offset[cabinet_idx][2]]
                        self.mesh = Cuboid(self.cab_up_down_inner_sizes[cabinet_idx],
                                           self.cabinet_size[cabinet_idx][0] - 2 * self.cab_left_right_inner_sizes[cabinet_idx],
                                           self.cabinet_size[cabinet_idx][2],
                                           position=mesh_position)
                    else:
                        mesh_position = [
                            self.cabinet_offset[cabinet_idx][0],
                            self.cabinet_offset[cabinet_idx][1] - self.cabinet_size[cabinet_idx][1] / 2,
                            self.cabinet_offset[cabinet_idx][2] - (self.cabinet_size[cabinet_idx][2] + self.cab_backs_size[cabinet_idx]) / 2]
                        self.mesh = Cuboid(self.cabinet_size[cabinet_idx][1],
                                           self.cabinet_size[cabinet_idx][0],
                                           self.cab_backs_size[cabinet_idx],
                                           position=mesh_position)
                else:
                    mesh_position = [
                        self.cabinet_offset[cabinet_idx][0],
                        - (self.layers_offset[cabinet_idx] + (mesh_idx - sum(total_mesh_list[:1])) * self.interval_between_layers[cabinet_idx]) - self.cabinet_size[cabinet_idx][1] / 2,
                        self.cabinet_offset[cabinet_idx][2]]
                    self.mesh = Cuboid(self.layers_sizes[cabinet_idx],
                                       self.cabinet_size[cabinet_idx][0] - 2 * self.cab_left_right_inner_sizes[cabinet_idx],
                                       self.cabinet_size[cabinet_idx][2],
                                       position=mesh_position)

                # special case for top and beneath cabinet to adjust the position of the cabinet 

                # Y-axis adjustment
                self.mesh.vertices[:, 1] += self.cabinet_size[cabinet_idx][1] / 2

                # X-axis adjustment
                if self.number_of_top_cabinet == 2:
                    if cabinet_idx % 2 == 0:
                        self.mesh.vertices[:, 0] += self.cabinet_size[cabinet_idx][0] / 2
                    else:
                        self.mesh.vertices[:, 0] -= self.cabinet_size[cabinet_idx][0] / 2

                vertices_list.append(self.mesh.vertices)
                faces_list.append(self.mesh.faces + total_num_vertices)
                total_num_vertices += len(self.mesh.vertices)

            # build doors and drawers
            for space_idx in range(self.number_of_layers[cabinet_idx] + 1):
                if self.number_of_layers[cabinet_idx] == 0:
                    _height = self.cabinet_size[cabinet_idx][1] - 2 * self.cab_up_down_inner_sizes[cabinet_idx]
                    _pos = 0
                elif space_idx == 0:
                    _height = self.layers_offset[actual_idx] - self.cab_up_down_inner_sizes[cabinet_idx] - self.layers_sizes[cabinet_idx] / 2
                    _pos = self.cabinet_size[cabinet_idx][1] - self.layers_offset[actual_idx] / 2 - self.cab_up_down_inner_sizes[cabinet_idx] / 2 + self.layers_sizes[cabinet_idx] / 4
                elif space_idx == self.number_of_layers[cabinet_idx]:
                    _height = self.cabinet_size[cabinet_idx][1] - (self.layers_offset[1] + (space_idx - 1) * self.interval_between_layers[1]) - \
                              self.cab_up_down_inner_sizes[cabinet_idx] - self.layers_sizes[cabinet_idx] / 2
                    _pos = self.cab_up_down_inner_sizes[cabinet_idx] + _height / 2 - 5
                else:
                    _height = self.interval_between_layers[0] - self.layers_sizes[cabinet_idx]
                    _pos = self.cabinet_size[cabinet_idx][1] - self.layers_offset[actual_idx] - (2 * space_idx - 1) / 2 * self.interval_between_layers[0]

                if self.type_of_spaces[cabinet_idx][space_idx] == 1:
                    for mesh_idx in range(5 + self.drawer_number_of_handles[cabinet_idx][space_idx]):
                        if mesh_idx < 2:
                            position_sign = -1 if mesh_idx == 0 else 1
                            mesh_position = [position_sign * (
                                    self.cabinet_size[cabinet_idx][0] / 2 - self.cab_left_right_inner_sizes[cabinet_idx] - self.drawer_inner_sizes[0] / 2) +
                                             self.cabinet_offset[cabinet_idx][0],
                                             _pos + self.cabinet_offset[cabinet_idx][1] - self.cabinet_size[cabinet_idx][1] / 2,
                                             self.cabinet_offset[cabinet_idx][2] + self.drawer_offset[cabinet_idx][space_idx] + self.drawer_interval[cabinet_idx][space_idx] / 2]
                            self.mesh = Cuboid(_height,
                                               self.drawer_inner_sizes[0],
                                               self.cabinet_size[cabinet_idx][2] - self.drawer_interval[cabinet_idx][space_idx],
                                               position=mesh_position)
                        elif mesh_idx < 4:
                            position_sign = -1 if mesh_idx == 3 else 1
                            mesh_position = [self.cabinet_offset[cabinet_idx][0],
                                             _pos + self.cabinet_offset[cabinet_idx][1] - self.cabinet_size[cabinet_idx][1] / 2,
                                             self.cabinet_offset[cabinet_idx][2] + position_sign * (self.cabinet_size[cabinet_idx][2] - self.drawer_inner_sizes[1]) / 2]
                            self.mesh = Cuboid(_height,
                                               self.cabinet_size[cabinet_idx][0] - 2 * self.cab_left_right_inner_sizes[cabinet_idx] - 2 * self.drawer_inner_sizes[0],
                                               self.drawer_inner_sizes[1],
                                               position=mesh_position)
                        elif mesh_idx == 4:
                            mesh_position = [self.cabinet_offset[cabinet_idx][0],
                                             _pos + self.cabinet_offset[cabinet_idx][1] - (_height + self.drawer_bottom_size[0]) / 2 - self.cabinet_size[cabinet_idx][1] / 2,
                                             self.cabinet_offset[cabinet_idx][2] + self.drawer_interval[cabinet_idx][space_idx] / 2 + self.drawer_offset[cabinet_idx][space_idx]]
                            self.mesh = Cuboid(self.drawer_bottom_size[0],
                                               self.cabinet_size[cabinet_idx][0] - 2 * self.cab_left_right_inner_sizes[cabinet_idx],
                                               self.cabinet_size[cabinet_idx][2] - self.drawer_interval[cabinet_idx][space_idx],
                                               position=mesh_position)
                        else:
                            if self.drawer_number_of_handles[cabinet_idx][space_idx] == 2:
                                position_sign = 1 if mesh_idx == 5 else -1
                            else:
                                position_sign = 0
                            mesh_rotation = [0, 0, self.drawer_handles_rotation[cabinet_idx][space_idx]]
                            mesh_position = [self.cabinet_offset[cabinet_idx][0] + self.drawer_handles_offsets[cabinet_idx][space_idx][0] +
                                             position_sign * self.drawer_handles_separation[cabinet_idx][space_idx],
                                             _pos + self.cabinet_offset[cabinet_idx][1] + self.drawer_handles_offsets[cabinet_idx][space_idx][1] - self.cabinet_size[cabinet_idx][1] / 2,
                                             self.cabinet_offset[cabinet_idx][2] + (self.cabinet_size[cabinet_idx][2] + self.drawer_handles_size[cabinet_idx][2]) / 2 +
                                             self.drawer_offset[cabinet_idx][space_idx]]
                            self.mesh = Cuboid(self.drawer_handles_size[cabinet_idx][1], self.drawer_handles_size[cabinet_idx][0], self.drawer_handles_size[cabinet_idx][2],
                                               position=mesh_position, rotation=mesh_rotation)

                        # special case for top and beneath cabinet to adjust the position of the cabinet 

                        # Y-axis adjustment
                        self.mesh.vertices[:, 1] += self.cabinet_size[cabinet_idx][1] / 2

                        # X-axis adjustment
                        if self.number_of_top_cabinet == 2:
                            if cabinet_idx % 2 == 0:
                                self.mesh.vertices[:, 0] += self.cabinet_size[cabinet_idx][0] / 2
                            else:
                                self.mesh.vertices[:, 0] -= self.cabinet_size[cabinet_idx][0] / 2

                        vertices_list.append(self.mesh.vertices)
                        faces_list.append(self.mesh.faces + total_num_vertices)
                        total_num_vertices += len(self.mesh.vertices)
                elif self.type_of_spaces[cabinet_idx][space_idx] == 2:
                    for mesh_idx in range(2):
                        if mesh_idx == 0:
                            mesh_position = [self.cabinet_offset[cabinet_idx][0],
                                             _pos + self.cabinet_offset[cabinet_idx][1],
                                             self.cabinet_offset[cabinet_idx][2] + (self.cabinet_size[cabinet_idx][2] + self.door_sizes[0]) / 2]
                            self.mesh = Cuboid(_height,
                                               self.cabinet_size[cabinet_idx][0] - 2 * self.cab_left_right_inner_sizes[cabinet_idx],
                                               self.door_sizes[0],
                                               position=mesh_position)
                        else:
                            mesh_rotation = [0, 0, self.door_handles_rotation[cabinet_idx][space_idx]]
                            mesh_position = [self.cabinet_offset[cabinet_idx][0] + self.door_handles_offsets[cabinet_idx][space_idx][0],
                                             _pos + self.cabinet_offset[cabinet_idx][1] + self.door_handles_offsets[cabinet_idx][space_idx][1],
                                             self.cabinet_offset[cabinet_idx][2] + (self.cabinet_size[cabinet_idx][2] + self.door_sizes[0]) / 2 + (
                                                     self.door_handles_size[cabinet_idx][2] + self.door_sizes[0]) / 2]
                            self.mesh = Cuboid(self.door_handles_size[1], self.door_handles_size[0], self.door_handles_size[2],
                                               position=mesh_position, rotation=mesh_rotation)
                    # special case for top and beneath cabinet to adjust the position of the cabinet 

                    # Y-axis adjustment
                    self.mesh.vertices[:, 1] += self.cabinet_size[cabinet_idx][1] / 2

                    # X-axis adjustment
                    if self.number_of_top_cabinet == 2:
                        if cabinet_idx % 2 == 0:
                            self.mesh.vertices[:, 0] += self.cabinet_size[cabinet_idx][0] / 2
                        else:
                            self.mesh.vertices[:, 0] -= self.cabinet_size[cabinet_idx][0] / 2
                        vertices_list.append(self.mesh.vertices)
                        faces_list.append(self.mesh.faces + total_num_vertices)
                        total_num_vertices += len(self.mesh.vertices)
                else:
                    pass

        # beneath cabinet
        for cabinet_idx in range(self.number_of_beneath_cabinet):
            actual_idx = 2 + cabinet_idx
            cabinet_body_mesh = 5
            layers_mesh = number_of_layers[actual_idx]
            cabinet_doors_mesh = sum(2 for type_of_space in self.type_of_spaces[cabinet_idx] if type_of_space == 2)
            cabinet_drawers_mesh = sum(
                (5 + self.drawer_number_of_handles[cabinet_idx][space_idx]) for space_idx, type_of_space in enumerate(self.type_of_spaces[cabinet_idx]) if type_of_space == 1)
            total_mesh_list = [cabinet_body_mesh, layers_mesh, cabinet_doors_mesh, cabinet_drawers_mesh]

            # build cabinet body and layers
            for mesh_idx in range(sum(total_mesh_list[:2])):
                if mesh_idx < sum(total_mesh_list[:1]):
                    position_sign = -1 if mesh_idx % 2 == 0 else 1
                    if mesh_idx < 2:
                        mesh_position = [
                            position_sign * (self.cabinet_size[cabinet_idx][0] - self.cab_left_right_inner_sizes[actual_idx]) / 2 + self.cabinet_offset[cabinet_idx][0],
                            self.cabinet_offset[cabinet_idx][1],
                            self.cabinet_offset[cabinet_idx][2]]
                        self.mesh = Cuboid(self.cabinet_size[cabinet_idx][1], self.cab_left_right_inner_sizes[actual_idx], self.cabinet_size[cabinet_idx][2],
                                           position=mesh_position)
                    elif mesh_idx < 4:
                        mesh_position = [
                            self.cabinet_offset[cabinet_idx][0],
                            position_sign * (self.cabinet_size[cabinet_idx][1] - self.cab_up_down_inner_sizes[actual_idx]) / 2 + self.cabinet_offset[cabinet_idx][1],
                            self.cabinet_offset[cabinet_idx][2]]
                        self.mesh = Cuboid(self.cab_up_down_inner_sizes[actual_idx],
                                           self.cabinet_size[cabinet_idx][0] - 2 * self.cab_left_right_inner_sizes[actual_idx],
                                           self.cabinet_size[cabinet_idx][2],
                                           position=mesh_position)
                    else:
                        if self.cab_backs_size[actual_idx] != 0:
                            mesh_position = [
                                self.cabinet_offset[cabinet_idx][0],
                                self.cabinet_offset[cabinet_idx][1],
                                self.cabinet_offset[cabinet_idx][2] - (self.cabinet_size[cabinet_idx][2] + self.cab_backs_size[actual_idx]) / 2]
                            self.mesh = Cuboid(self.cabinet_size[cabinet_idx][1],
                                               self.cabinet_size[cabinet_idx][0],
                                               self.cab_backs_size[actual_idx],
                                               position=mesh_position)
                        else:
                            pass
                else:
                    mesh_position = [
                        self.cabinet_offset[cabinet_idx][0],
                        self.cabinet_offset[cabinet_idx][1] + self.cabinet_size[cabinet_idx][1] / 2 - (
                                self.layers_offset[actual_idx] + (mesh_idx - sum(total_mesh_list[:1])) * self.interval_between_layers[actual_idx]),
                        self.cabinet_offset[cabinet_idx][2]]
                    self.mesh = Cuboid(self.layers_sizes[actual_idx],
                                       self.cabinet_size[cabinet_idx][0] - 2 * self.cab_left_right_inner_sizes[actual_idx],
                                       self.cabinet_size[cabinet_idx][2],
                                       position=mesh_position)

                vertices_list.append(self.mesh.vertices)
                faces_list.append(self.mesh.faces + total_num_vertices)
                total_num_vertices += len(self.mesh.vertices)

            # build doors and drawers
            for space_idx in range(self.number_of_layers[actual_idx] + 1):
                if self.number_of_layers[actual_idx] == 0:
                    _height = self.cabinet_size[cabinet_idx][1] - 2 * self.cab_up_down_inner_sizes[actual_idx]
                    _pos = self.cabinet_size[cabinet_idx][1] / 2
                elif space_idx == 0:
                    _height = self.layers_offset[actual_idx] - self.cab_up_down_inner_sizes[actual_idx] - self.layers_sizes[actual_idx] / 2
                    _pos = self.cabinet_size[cabinet_idx][1] - self.layers_offset[actual_idx] / 2 - self.cab_up_down_inner_sizes[actual_idx] / 2 + self.layers_sizes[actual_idx] / 4
                elif space_idx == self.number_of_layers[actual_idx]:
                    _height = self.cabinet_size[cabinet_idx][1] - (self.layers_offset[actual_idx] + (space_idx - 1) * self.interval_between_layers[0]) - \
                              self.cab_up_down_inner_sizes[actual_idx] - self.layers_sizes[actual_idx] / 2
                    _pos = self.cab_up_down_inner_sizes[actual_idx] + _height / 2
                else:
                    _height = self.interval_between_layers[0] - self.layers_sizes[actual_idx]
                    _pos = self.cabinet_size[cabinet_idx][1] - self.layers_offset[actual_idx] - (2 * space_idx - 1) / 2 * self.interval_between_layers[0]

                if self.type_of_spaces[cabinet_idx][space_idx] == 1:
                    for mesh_idx in range(5 + self.drawer_number_of_handles[cabinet_idx][space_idx]):
                        if mesh_idx < 2:
                            position_sign = -1 if mesh_idx == 0 else 1
                            mesh_position = [position_sign * (
                                    self.cabinet_size[cabinet_idx][0] / 2 - self.cab_left_right_inner_sizes[actual_idx] - self.drawer_inner_sizes[0] / 2) +
                                             self.cabinet_offset[cabinet_idx][0],
                                             _pos + self.cabinet_offset[cabinet_idx][1],
                                             self.cabinet_offset[cabinet_idx][2] + self.drawer_offset[cabinet_idx][space_idx] + self.drawer_interval[cabinet_idx][space_idx] / 2]
                            self.mesh = Cuboid(_height,
                                               self.drawer_inner_sizes[0],
                                               self.cabinet_size[cabinet_idx][2] - self.drawer_interval[cabinet_idx][space_idx],
                                               position=mesh_position)
                        elif mesh_idx < 4:
                            position_sign = -1 if mesh_idx == 3 else 1
                            mesh_position = [self.cabinet_offset[cabinet_idx][0],
                                             _pos + self.cabinet_offset[cabinet_idx][1],
                                             self.cabinet_offset[cabinet_idx][2] + position_sign * (self.cabinet_size[cabinet_idx][2] - self.drawer_inner_sizes[1]) / 2 + self.drawer_offset[cabinet_idx][space_idx]]
                            self.mesh = Cuboid(_height,
                                               self.cabinet_size[cabinet_idx][0] - 2 * self.cab_left_right_inner_sizes[actual_idx] - 2 * self.drawer_inner_sizes[0],
                                               self.drawer_inner_sizes[1],
                                               position=mesh_position)
                        elif mesh_idx == 4:
                            mesh_position = [self.cabinet_offset[cabinet_idx][0],
                                             _pos + self.cabinet_offset[cabinet_idx][1] - (_height + self.drawer_bottom_size[0]) / 2,
                                             self.cabinet_offset[cabinet_idx][2] + self.drawer_interval[cabinet_idx][space_idx] / 2 + self.drawer_offset[cabinet_idx][space_idx]]
                            self.mesh = Cuboid(self.drawer_bottom_size[0],
                                               self.cabinet_size[cabinet_idx][0] - 2 * self.cab_left_right_inner_sizes[actual_idx],
                                               self.cabinet_size[cabinet_idx][2] - self.drawer_interval[cabinet_idx][space_idx],
                                               position=mesh_position)
                        else:
                            if self.drawer_number_of_handles[cabinet_idx][space_idx] == 2:
                                position_sign = 1 if mesh_idx == 5 else -1
                            else:
                                position_sign = 0
                            mesh_rotation = [0, 0, self.drawer_handles_rotation[cabinet_idx][space_idx]]
                            mesh_position = [self.cabinet_offset[cabinet_idx][0] + self.drawer_handles_offsets[cabinet_idx][space_idx][0] +
                                             position_sign * self.drawer_handles_separation[cabinet_idx][space_idx],
                                             _pos + self.cabinet_offset[cabinet_idx][1] + self.drawer_handles_offsets[cabinet_idx][space_idx][1],
                                             self.cabinet_offset[cabinet_idx][2] + (self.cabinet_size[cabinet_idx][2] + self.drawer_handles_size[cabinet_idx][2]) / 2 +
                                             self.drawer_offset[cabinet_idx][space_idx]]
                            self.mesh = Cuboid(self.drawer_handles_size[cabinet_idx][1], self.drawer_handles_size[cabinet_idx][0], self.drawer_handles_size[cabinet_idx][2],
                                               position=mesh_position, rotation=mesh_rotation)

                        # special case for top and beneath cabinet to adjust the position of the cabinet 

                        # Y-axis adjustment
                        self.mesh.vertices[:, 1] -= self.cabinet_size[cabinet_idx][1] / 2

                        vertices_list.append(self.mesh.vertices)
                        faces_list.append(self.mesh.faces + total_num_vertices)
                        total_num_vertices += len(self.mesh.vertices)
                elif self.type_of_spaces[cabinet_idx][space_idx] == 2:
                    for mesh_idx in range(2):
                        if mesh_idx == 0:
                            mesh_position = [self.cabinet_offset[cabinet_idx][0],
                                             _pos + self.cabinet_offset[cabinet_idx][1],
                                             self.cabinet_offset[cabinet_idx][2] + (self.cabinet_size[cabinet_idx][2] + self.door_sizes[0]) / 2]
                            self.mesh = Cuboid(_height,
                                               self.cabinet_size[cabinet_idx][0] - 2 * self.cab_left_right_inner_sizes[actual_idx],
                                               self.door_sizes[0],
                                               position=mesh_position)
                        else:
                            mesh_rotation = [0, 0, self.door_handles_rotation[cabinet_idx][space_idx]]
                            mesh_position = [self.cabinet_offset[cabinet_idx][0] + self.door_handles_offsets[cabinet_idx][space_idx][0],
                                             _pos + self.cabinet_offset[cabinet_idx][1] + self.door_handles_offsets[cabinet_idx][space_idx][1],
                                             self.cabinet_offset[cabinet_idx][2] + (self.cabinet_size[cabinet_idx][2] + self.door_sizes[0]) / 2 + (
                                                     self.door_handles_size[cabinet_idx][2] + self.door_sizes[0]) / 2]
                            self.mesh = Cuboid(self.door_handles_size[cabinet_idx][1], self.door_handles_size[cabinet_idx][0], self.door_handles_size[cabinet_idx][2],
                                               position=mesh_position, rotation=mesh_rotation)

                    # special case for top and beneath cabinet to adjust the position of the cabinet 

                        # Y-axis adjustment
                        self.mesh.vertices[:, 1] -= self.cabinet_size[cabinet_idx][1] / 2

                        vertices_list.append(self.mesh.vertices)
                        faces_list.append(self.mesh.faces + total_num_vertices)
                        total_num_vertices += len(self.mesh.vertices)
                else:
                    pass

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Cabinet'


class Regular_partition(ConceptTemplate):
    def __init__(self, has_partition, left_right_size, rear_size, left_right_separation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.has_partition = has_partition
        self.left_right_size = left_right_size
        self.rear_size = rear_size
        self.left_right_separation = left_right_separation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        if has_partition[0] == 1:
            mesh_position = [left_right_separation[0] / 2 + left_right_size[0] / 2, left_right_size[1] / 2, 0]
            self.mesh = Cuboid(left_right_size[1], left_right_size[0], left_right_size[2],
                               position=mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        if has_partition[1] == 1:
            mesh_position = [0, rear_size[0] / 2, -left_right_size[2] / 2 - rear_size[1] / 2]
            self.mesh = Cuboid(rear_size[0], left_right_separation[0] + 2 * left_right_size[0], rear_size[1],
                               position=mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        if has_partition[2] == 1:
            mesh_position = [-left_right_separation[0] / 2 - left_right_size[0] / 2, left_right_size[1] / 2, 0]
            self.mesh = Cuboid(left_right_size[1], left_right_size[0], left_right_size[2],
                               position=mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Board'