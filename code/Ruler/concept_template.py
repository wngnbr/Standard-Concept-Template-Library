import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh


class Symmetrical_body(ConceptTemplate):
    def __init__(self, size, separation, left_right_offset, body_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        body_rotation = [x / 180 * np.pi for x in body_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.separation = separation  # offset[0]
        self.left_right_offset = left_right_offset
        self.body_rotation = body_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # left, more rotation more down
        mesh_1_rotation = [0, 0, body_rotation[0]]
        mesh_1_position = [-separation[0] - (size[0] * np.cos(body_rotation[0]) + size[1] * np.sin(body_rotation[0])) / 2,
                           (size[1] * np.cos(body_rotation[0]) - size[0] * np.sin(body_rotation[0])) / 2, -left_right_offset[0]]
        self.mesh_1 = Cuboid(size[1], size[0], size[2], position=mesh_1_position, rotation=mesh_1_rotation)
        vertices_list.append(self.mesh_1.vertices)
        faces_list.append(self.mesh_1.faces + total_num_vertices)
        total_num_vertices += len(self.mesh_1.vertices)

        mesh_2_rotation = [0, 0, -body_rotation[0]]
        mesh_2_position = [separation[0] + (size[0] * np.cos(body_rotation[0]) + size[1] * np.sin(body_rotation[0])) / 2,
                           (size[1] * np.cos(body_rotation[0]) - size[0] * np.sin(body_rotation[0])) / 2, left_right_offset[0]]
        self.mesh_2 = Cuboid(size[1], size[0], size[2], position=mesh_2_position, rotation=mesh_2_rotation)
        vertices_list.append(self.mesh_2.vertices)
        faces_list.append(self.mesh_2.faces + total_num_vertices)
        total_num_vertices += len(self.mesh_2.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Ruler'


class Asymmetrical_body(ConceptTemplate):
    def __init__(self, left_size, right_size, separation, left_right_offset, body_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        body_rotation = [x / 180 * np.pi for x in body_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.left_size = left_size
        self.right_size = right_size
        self.separation = separation  # offset[0]
        self.left_right_offset = left_right_offset
        self.body_rotation = body_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_1_rotation = [0, 0, body_rotation[0]]
        mesh_1_position = [-separation[0] - left_size[0] / 2 * np.cos(body_rotation[0]) - left_size[1] / 2 * np.sin(body_rotation[0]),
                           -left_size[0] / 2 * np.sin(body_rotation[0]) + left_size[1] / 2 * np.cos(body_rotation[0]), -left_right_offset[0]]
        self.mesh_1 = Cuboid(left_size[1], left_size[0], left_size[2], position=mesh_1_position, rotation=mesh_1_rotation)
        vertices_list.append(self.mesh_1.vertices)
        faces_list.append(self.mesh_1.faces + total_num_vertices)
        total_num_vertices += len(self.mesh_1.vertices)

        mesh_2_rotation = [0, 0, -body_rotation[0]]
        mesh_2_position = [separation[0] + right_size[0] / 2 * np.cos(body_rotation[0]) - right_size[1] / 2 * np.sin(-body_rotation[0]),
                           right_size[0] / 2 * np.sin(-body_rotation[0]) + right_size[1] / 2 * np.cos(body_rotation[0]), left_right_offset[0]]
        self.mesh_2 = Cuboid(right_size[1], right_size[0], right_size[2], position=mesh_2_position, rotation=mesh_2_rotation)
        vertices_list.append(self.mesh_2.vertices)
        faces_list.append(self.mesh_2.faces + total_num_vertices)
        total_num_vertices += len(self.mesh_2.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Ruler'


class Regular_shaft(ConceptTemplate):
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

        mesh_rotation = [np.pi / 2, 0, 0]
        self.mesh = Cylinder(size[1], size[0], rotation=mesh_rotation)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Shaft'