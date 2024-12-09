import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh

class Standard_Body(ConceptTemplate):
    def __init__(self, base_size, beside_size, beside_seperation, beside_offset_z, has_shaft, shaft_central_size, shaft_beside_size, shaft_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.base_size = base_size
        self.beside_size = beside_size
        self.beside_seperation = beside_seperation
        self.beside_offset_z = beside_offset_z
        self.has_shaft = has_shaft
        self.shaft_central_size = shaft_central_size
        self.shaft_beside_size = shaft_beside_size
        self.shaft_offset = shaft_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.base_mesh = Cuboid(base_size[1], base_size[0], base_size[2])
        vertices_list.append(self.base_mesh.vertices)
        faces_list.append(self.base_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.base_mesh.vertices)

        left_mesh_position = [
            beside_seperation[0] / 2,
            (base_size[1] + beside_size[1]) / 2,
            -base_size[2] / 2 + beside_size[2] / 2 + beside_offset_z[0]
        ]
        self.left_mesh = Cuboid(beside_size[1], beside_size[0], beside_size[2],
                                position = left_mesh_position)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            -beside_seperation[0] / 2,
            (base_size[1] + beside_size[1]) / 2,
            -base_size[2] / 2 + beside_size[2] / 2 + beside_offset_z[0]
        ]
        self.right_mesh = Cuboid(beside_size[1], beside_size[0], beside_size[2],
                                 position = right_mesh_position)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        if has_shaft[0] == 1:
            central_shaft_mesh_position = [
                0,
                (base_size[1] + beside_size[1]) / 2 + shaft_offset[0],
                -base_size[2] / 2 + beside_size[2] / 2 + beside_offset_z[0] + shaft_offset[1]
            ]
            central_shaft_mesh_rotation = [0, 0, np.pi / 2]
            self.central_shaft_mesh = Cylinder(shaft_central_size[1], shaft_central_size[0],
                                               position = central_shaft_mesh_position,
                                               rotation = central_shaft_mesh_rotation)
            vertices_list.append(self.central_shaft_mesh.vertices)
            faces_list.append(self.central_shaft_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.central_shaft_mesh.vertices)

            left_shaft_mesh_position = [
                (shaft_central_size[1] + shaft_beside_size[1]) / 2,
                (base_size[1] + beside_size[1]) / 2 + shaft_offset[0],
                -base_size[2] / 2 + beside_size[2] / 2 + beside_offset_z[0] + shaft_offset[1]
            ]
            left_shaft_mesh_rotation = [0, 0, np.pi / 2]
            self.left_shaft_mesh = Cylinder(shaft_beside_size[1], shaft_beside_size[0],
                                            position = left_shaft_mesh_position,
                                            rotation = left_shaft_mesh_rotation)
            vertices_list.append(self.left_shaft_mesh.vertices)
            faces_list.append(self.left_shaft_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.left_shaft_mesh.vertices)

            right_shaft_mesh_position = [
                -(shaft_central_size[1] + shaft_beside_size[1]) / 2,
                (base_size[1] + beside_size[1]) / 2 + shaft_offset[0],
                -base_size[2] / 2 + beside_size[2] / 2 + beside_offset_z[0] + shaft_offset[1]
            ]
            right_shaft_mesh_rotation = [0, 0, np.pi / 2]
            self.right_shaft_mesh = Cylinder(shaft_beside_size[1], shaft_beside_size[0],
                                             position = right_shaft_mesh_position,
                                             rotation = right_shaft_mesh_rotation)
            vertices_list.append(self.right_shaft_mesh.vertices)
            faces_list.append(self.right_shaft_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.right_shaft_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Simplified_Cover(ConceptTemplate):
    def __init__(self, size, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0,
            size[1] / 2,
            size[2] / 2
        ]
        self.mesh = Cuboid(size[1], size[0], size[2],
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

        self.semantic = 'Cover'


class Carved_Cover(ConceptTemplate):
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


class Carved_Magazine(ConceptTemplate):
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
        top_mesh_rotation = [0, 0, np.pi]
        self.top_mesh = Cuboid(outer_size[1] - inner_size[1], outer_size[0], outer_size[2],
                               position = top_mesh_position,
                               rotation = top_mesh_rotation)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0, 
            inner_size[1] / 2,
            outer_size[2] / 2
        ]
        bottom_mesh_rotation = [0, 0, np.pi]
        self.bottom_mesh = Rectangular_Ring(inner_size[1], outer_size[0], outer_size[2], inner_size[0], inner_size[2],
                                            position = bottom_mesh_position,
                                            rotation = bottom_mesh_rotation)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Magazine'


class Complex_Magazine(ConceptTemplate):
    def __init__(self, size, thickness, front_height, beside_length, beside_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.thickness = thickness
        self.front_height = front_height
        self.beside_length = beside_length
        self.beside_offset = beside_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            0,
            thickness[0] / 2,
            size[2] / 2
        ]
        self.bottom_mesh = Cuboid(thickness[0], size[0], size[2],
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        behind_mesh_position = [
            0,
            size[1] / 2,
            thickness[0] / 2
        ]
        self.behind_mesh = Cuboid(size[1], size[0], thickness[0],
                                  position = bottom_mesh_position)
        vertices_list.append(self.behind_mesh.vertices)
        faces_list.append(self.behind_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.behind_mesh.vertices)

        front_mesh_position = [
            0,
            front_height[0] / 2,
            size[2] - thickness[0] / 2
        ]
        self.front_mesh = Cuboid(front_height[0], size[0], thickness[0],
                                  position = front_mesh_position)
        vertices_list.append(self.front_mesh.vertices)
        faces_list.append(self.front_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_mesh.vertices)

        left_mesh_position = [
            size[0] / 2 - thickness[0] / 2,
            size[1] / 2,
            size[2] / 2 + beside_offset[0]
        ]
        self.left_mesh = Cuboid(size[1], thickness[0], beside_length[0],
                                position = left_mesh_position)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            -size[0] / 2 + thickness[0] / 2,
            size[1] / 2,
            size[2] / 2 + beside_offset[0]
        ]
        self.right_mesh = Cuboid(size[1], thickness[0], beside_length[0],
                                position = right_mesh_position)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Magazine'