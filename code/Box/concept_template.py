import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh

class Cuboidal_Body(ConceptTemplate):
    def __init__(self, top_size, bottom_size, height, top_bottom_offset, thickness, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.top_size = top_size
        self.bottom_size = bottom_size
        self.height = height
        self.top_bottom_offset = top_bottom_offset
        self.thickness = thickness

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        middle_x = top_size[0] * thickness[1] / height[0] + bottom_size[0] * (height[0] - thickness[1]) / height[0]
        middle_z = top_size[1] * thickness[1] / height[0] + bottom_size[1] * (height[0] - thickness[1]) / height[0]
        middle_offset_x = top_bottom_offset[0] * thickness[1] / height[0]
        middle_offset_z = top_bottom_offset[1] * thickness[1] / height[0]


        top_mesh_position = [
            middle_offset_x,
            thickness[1] / 2,
            middle_offset_z
        ]
        self.top_mesh = Rectangular_Ring(height[0] - thickness[1], top_size[0], top_size[1],
                                         top_size[0] - thickness[0] * 2, top_size[1] - thickness[2] * 2,
                                         [0, 0], middle_x, middle_z,
                                         middle_x - thickness[0] * 2, middle_z - thickness[2] * 2,
                                         top_bottom_offset = [top_bottom_offset[0] - middle_offset_x, top_bottom_offset[1] - middle_offset_z],
                                         position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0, 
            -(height[0] - thickness[1]) / 2,
            0
        ]
        self.bottom_mesh = Cuboid(thickness[1], middle_x, middle_z,
                                  bottom_size[0], bottom_size[1],
                                  top_offset = [middle_offset_x, middle_offset_z],
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)


        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Fourfold_Cover(ConceptTemplate):
    def __init__(self, has_cover, front_behind_size, left_right_size, cover_separation, cover_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        cover_rotation = [x / 180 * np.pi for x in cover_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.has_cover = has_cover
        self.front_behind_size = front_behind_size
        self.left_right_size = left_right_size
        self.cover_separation = cover_separation
        self.cover_rotation = cover_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        if has_cover[0] == 1:
            mesh_position = [
                0, 
                front_behind_size[1] * np.cos(cover_rotation[0]) / 2, 
                cover_separation[0] / 2 + front_behind_size[1] * np.sin(cover_rotation[0]) / 2
            ]
            mesh_rotation = [
                cover_rotation[0], 
                0, 
                0
            ]
            self.mesh = Cuboid(front_behind_size[1], front_behind_size[0], front_behind_size[2],
                               position = mesh_position,
                               rotation = mesh_rotation)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        if has_cover[1] == 1:
            mesh_position = [
                0, 
                front_behind_size[1] * np.cos(cover_rotation[1]) / 2, 
                -cover_separation[0] / 2 + front_behind_size[1] * np.sin(cover_rotation[1]) / 2
            ]
            mesh_rotation = [
                cover_rotation[1], 
                0, 
                0
            ]
            self.mesh = Cuboid(front_behind_size[1], front_behind_size[0], front_behind_size[2],
                               position = mesh_position,
                               rotation = mesh_rotation)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        if has_cover[2] == 1:
            mesh_position = [
                -cover_separation[1] / 2 - left_right_size[1] * np.sin(cover_rotation[2]) / 2, 
                left_right_size[1] * np.cos(cover_rotation[2]) / 2, 
                0
            ]
            mesh_rotation = [
                0, 
                0, 
                cover_rotation[2]
            ]
            self.mesh = Cuboid(left_right_size[1], left_right_size[0], left_right_size[2],
                               position = mesh_position,
                               rotation = mesh_rotation)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        if has_cover[3] == 1:
            mesh_position = [
                cover_separation[1] / 2 - left_right_size[1] * np.sin(cover_rotation[3]) / 2, 
                left_right_size[1] * np.cos(cover_rotation[3]) / 2, 
                0
            ]
            mesh_rotation = [
                0, 
                0, 
                cover_rotation[3]
            ]
            self.mesh = Cuboid(left_right_size[1], left_right_size[0], left_right_size[2],
                               position = mesh_position,
                               rotation = mesh_rotation)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)


        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Cover'


class Regular_Cover(ConceptTemplate):
    def __init__(self, outer_size, inner_size, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.outer_size = outer_size
        self.inner_size = inner_size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        top_mesh_position = [
            0,
            (outer_size[1] + inner_size[1]) / 2,
            outer_size[2] / 2
        ]
        self.top_mesh = Cuboid(outer_size[1] - inner_size[1], outer_size[0], outer_size[2],
                               position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0, 
            inner_size[1] / 2,
            outer_size[2] / 2
        ]
        self.bottom_mesh = Rectangular_Ring(inner_size[1], outer_size[0], outer_size[2], inner_size[0], inner_size[2],
                                            position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)


        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Cover'

class Cuboidal_Leg(ConceptTemplate):
    def __init__(self, front_legs_size, rear_legs_size, legs_separation, num_legs, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.front_legs_size = front_legs_size
        self.rear_legs_size = rear_legs_size
        self.legs_separation = legs_separation
        self.num_legs = num_legs

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        if (num_legs[0] == 1):
            mesh_position = [
                0,
                -front_legs_size[1] / 2,
                0
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        elif (num_legs[0] == 2):
            mesh_position = [
                legs_separation[0] / 2, 
                -front_legs_size[1] / 2,
                0
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

            mesh_position = [
                -legs_separation[0] / 2,
                -front_legs_size[1] / 2,
                0
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        elif (num_legs[0] == 3):
            mesh_position = [
                legs_separation[0] / 2,
                -front_legs_size[1] / 2,
                legs_separation[2] / 2
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

            mesh_position = [
                -legs_separation[0] / 2,
                -front_legs_size[1] / 2,
                legs_separation[2] / 2
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

            mesh_position = [
                0,
                -rear_legs_size[1] / 2,
                -legs_separation[2] / 2
            ]
            self.mesh = Cuboid(rear_legs_size[1], rear_legs_size[0], rear_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        elif (num_legs[0] == 4):
            mesh_position = [
                legs_separation[0] / 2,
                -front_legs_size[1] / 2,
                legs_separation[2] / 2
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

            mesh_position = [
                -legs_separation[0] / 2,
                -front_legs_size[1] / 2,
                legs_separation[2] / 2
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

            mesh_position = [
                legs_separation[1] / 2,
                -rear_legs_size[1] / 2,
                -legs_separation[2] / 2
            ]
            self.mesh = Cuboid(rear_legs_size[1], rear_legs_size[0], rear_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

            mesh_position = [
                -legs_separation[1] / 2,
                -rear_legs_size[1] / 2,
                -legs_separation[2] / 2
            ]
            self.mesh = Cuboid(rear_legs_size[1], rear_legs_size[0], rear_legs_size[2],
                               position = mesh_position)
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