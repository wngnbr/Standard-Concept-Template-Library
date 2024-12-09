import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh


class Regular_Body(ConceptTemplate):
    def __init__(self, bottom_size, middle_size, top_size, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.bottom_size = bottom_size
        self.middle_size = middle_size
        self.top_size = top_size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        base_mesh_position = [0, (top_size[0] - bottom_size[1]) / 2, 0]
        base_mesh_rotation = [0, 0, -np.pi / 2]
        self.base_mesh = Rectangular_Ring(
            bottom_size[0],
            top_size[0] + middle_size[0] + bottom_size[1],
            bottom_size[2],
            middle_size[0],
            bottom_size[2] - middle_size[1] * 2,
            inner_offset=[(top_size[0] - bottom_size[1]) / 2, 0],
            top_bottom_offset=[0, 0],
            position=base_mesh_position,
            rotation=base_mesh_rotation,
        )
        vertices_list.append(self.base_mesh.vertices)
        faces_list.append(self.base_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.base_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Regular_Blade(ConceptTemplate):
    def __init__(self, bottom_size, top_size, top_bottom_offset, root_offset, blade_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        blade_rotation = [-x / 180 * np.pi for x in blade_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.bottom_size = bottom_size
        self.top_size = top_size
        self.top_bottom_offset = top_bottom_offset
        self.root_offset = root_offset
        self.blade_rotation = blade_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            bottom_size[0] / 2 - root_offset[0], 
            -bottom_size[1] / 2, 
            0
        ]
        self.bottom_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                                  position=bottom_mesh_position)
        self.bottom_mesh.vertices = apply_transformation(self.bottom_mesh.vertices, [0, 0, 0], [0, 0, blade_rotation[0]])
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [
            bottom_size[0] + top_size[0] / 2 - root_offset[0],
            top_bottom_offset[0] - bottom_size[1] / 2,
            0,
        ]
        self.top_mesh = Cuboid(top_size[1], top_size[0], top_size[2],
                               position=top_mesh_position)
        self.top_mesh.vertices = apply_transformation(self.top_mesh.vertices, [0, 0, 0], [0, 0, blade_rotation[0]])
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Blade'