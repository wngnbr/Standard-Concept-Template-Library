import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh


class Cylindrical_Body(ConceptTemplate):
    def __init__(self, size, x_z_ratio, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.x_z_ratio = x_z_ratio

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.base_mesh = Cylinder(size[2], size[0], size[1],
                                  x_z_ratio[0] * size[0], x_z_ratio[0] * size[1])
        vertices_list.append(self.base_mesh.vertices)
        faces_list.append(self.base_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.base_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Toothpaste_Body(ConceptTemplate):
    def __init__(self, radius, bottom_length, height, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.radius = radius
        self.bottom_length = bottom_length
        self.height = height

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.base_mesh = Cylinder(height[0], radius[0],
                                  bottom_radius=0,
                                  bottom_radius_z=bottom_length[0] / 2)
        vertices_list.append(self.base_mesh.vertices)
        faces_list.append(self.base_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.base_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Cuboidal_Body(ConceptTemplate):
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

        self.base_mesh = Cuboid(size[2], size[0], size[1])
        vertices_list.append(self.base_mesh.vertices)
        faces_list.append(self.base_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.base_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Domed_Cover(ConceptTemplate):
    def __init__(self, bottom_size, sphere_radius, sphere_exist_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        sphere_exist_rotation = [x / 180 * np.pi for x in sphere_exist_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.bottom_size = bottom_size
        self.sphere_radius = sphere_radius
        self.sphere_exist_rotation = sphere_exist_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            0,
            bottom_size[2] / 2,
            0,
        ]
        self.bottom_mesh = Cylinder(bottom_size[2], bottom_size[0], bottom_size[1],
                                    position=bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [
            0,
            bottom_size[2] - sphere_radius[1] * np.cos(sphere_exist_rotation[0]),
            0,
        ]
        self.top_mesh = Sphere(
            radius=sphere_radius[0],
            top_angle=0,
            bottom_angle=sphere_exist_rotation[0],
            radius_y=sphere_radius[1],
            position=top_mesh_position,
        )
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Cover'


class Cylindrical_Cover(ConceptTemplate):
    def __init__(self, outer_size, inner_size, position=[0, 0, 0], rotation=[0, 0, 0]):

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

        radius_middle = outer_size[0] + (outer_size[2] - inner_size[2]) * (outer_size[1] - outer_size[0]) / outer_size[2]

        ring_mesh_position = [0, inner_size[2] / 2, 0]
        self.ring_mesh = Ring(
            height=inner_size[2],
            outer_top_radius=radius_middle,
            inner_top_radius=inner_size[0],
            outer_bottom_radius=outer_size[1],
            inner_bottom_radius=inner_size[1],
            position=ring_mesh_position,
        )
        vertices_list.append(self.ring_mesh.vertices)
        faces_list.append(self.ring_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.ring_mesh.vertices)

        top_mesh_position = [
            0, 
            (inner_size[2] + outer_size[2]) / 2, 
            0
        ]
        self.top_mesh = Cylinder(outer_size[2] - inner_size[2], outer_size[0], radius_middle,
                                 position=top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)
        
        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Cover'
