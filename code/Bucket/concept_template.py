import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh

class Cylindrical_Body(ConceptTemplate):
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

        middle_radius = outer_size[0] * (1 - inner_size[2] / outer_size[2]) + outer_size[1] * inner_size[2] / outer_size[2]
        bottom_height = outer_size[2] - inner_size[2]
        bottom_mesh_position = [0, -inner_size[2] / 2, 0]
        self.bottom_mesh = Cylinder(bottom_height, middle_radius, outer_size[1], 
                                    position=bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [0, (outer_size[2] - inner_size[2]) / 2, 0]
        self.top_mesh = Ring(inner_size[2], outer_size[0], inner_size[0], 
                             outer_bottom_radius = middle_radius,
                             inner_bottom_radius = inner_size[1],
                             position=top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)


        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Prismatic_Body(ConceptTemplate):
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

        middle_radius = outer_size[0] * (1 - inner_size[2] / outer_size[2]) + outer_size[1] * inner_size[2] / outer_size[2]
        bottom_height = outer_size[2] - inner_size[2]
        bottom_mesh_position = [0, -inner_size[2] / 2, 0]
        self.bottom_mesh = Cuboid(bottom_height, middle_radius * np.sqrt(2), middle_radius * np.sqrt(2), 
                                  outer_size[1] * np.sqrt(2), outer_size[1] * np.sqrt(2),
                                  position=bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [0, (outer_size[2] - inner_size[2]) / 2, 0]
        self.top_mesh = Rectangular_Ring(inner_size[2], outer_size[0] * np.sqrt(2), outer_size[0] * np.sqrt(2), 
                                         inner_size[0] * np.sqrt(2), inner_size[0] * np.sqrt(2),
                                         outer_bottom_length = middle_radius * np.sqrt(2), outer_bottom_width = middle_radius * np.sqrt(2), 
                                         inner_bottom_length = inner_size[1] * np.sqrt(2), inner_bottom_width = inner_size[1] * np.sqrt(2), 
                                         position=top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)


        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Trifold_Handle(ConceptTemplate):
    def __init__(self, vertical_thickness, vertical_length, horizontal_thickness, vertical_rotation, vertical_separation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        vertical_rotation = [x / 180 * np.pi for x in vertical_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.vertical_thickness = vertical_thickness
        self.vertical_length = vertical_length
        self.horizontal_thickness = horizontal_thickness
        self.vertical_rotation = vertical_rotation
        self.vertical_separation = vertical_separation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        top_mesh_position = [
            vertical_separation[0] / 2 - vertical_length[0] * np.sin(vertical_rotation[0]) / 2, 
            vertical_length[0] * np.cos(vertical_rotation[0]) / 2, 
            0
            ]
        top_mesh_rotation = [0, 0, vertical_rotation[0]]
        self.top_mesh = Cuboid(vertical_length[0], vertical_thickness[0], vertical_thickness[1], 
                               position=top_mesh_position,
                               rotation=top_mesh_rotation)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            -vertical_separation[0] / 2 - vertical_length[1] * np.sin(vertical_rotation[1]) / 2, 
            vertical_length[1] * np.cos(vertical_rotation[1]) / 2, 
            0
            ]
        bottom_mesh_rotation = [0, 0, vertical_rotation[1]]
        self.bottom_mesh = Cuboid(vertical_length[1], vertical_thickness[0], vertical_thickness[1], 
                                  position=bottom_mesh_position,
                                  rotation=bottom_mesh_rotation)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        delta_x = vertical_separation[0] - vertical_length[0] * np.sin(vertical_rotation[0]) + vertical_length[1] * np.sin(vertical_rotation[1])
        delta_y = vertical_length[0] * np.cos(vertical_rotation[0]) - vertical_length[1] * np.cos(vertical_rotation[1])
        horizontal_length = np.sqrt(delta_y * delta_y + delta_x * delta_x) + vertical_thickness[0]
        horizontal_rotation = np.arctan(delta_y / delta_x)
        vertical_x_offset = (-vertical_length[0] * np.sin(vertical_rotation[0]) - vertical_length[1] * np.sin(vertical_rotation[1])) / 2
        vertical_y_offset = (vertical_length[1] * np.cos(vertical_rotation[1]) + vertical_length[0] * np.cos(vertical_rotation[0])) / 2
        vertical_mesh_position = [
            vertical_x_offset, 
            vertical_y_offset + horizontal_thickness[0] / 2, 
            0
            ]
        vertical_mesh_rotation = [0, 0, horizontal_rotation]
        self.vertical_mesh = Cuboid(horizontal_thickness[0], horizontal_length, horizontal_thickness[1], 
                                    position=vertical_mesh_position,
                                    rotation=vertical_mesh_rotation)
        vertices_list.append(self.vertical_mesh.vertices)
        faces_list.append(self.vertical_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.vertical_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)


        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Curved_Handle(ConceptTemplate):
    def __init__(self, radius, exist_angle, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        exist_angle = [x / 180 * np.pi for x in exist_angle]
        super().__init__(position, rotation)

        # Record Parameters
        self.radius = radius
        self.exist_angle = exist_angle

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_rotation = [-np.pi / 2, 0, 0]
        self.mesh = Torus(radius[0], radius[1], exist_angle[0],
                          rotation = mesh_rotation)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order='ZXY')

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Round_U_Handle(ConceptTemplate):
    def __init__(self, inner_radius, vertical_separation, vertical_length, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.inner_radius = inner_radius
        self.vertical_separation = vertical_separation
        self.vertical_length = vertical_length

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        left_mesh_position = [
            vertical_separation[0] / 2, 
            vertical_length[0] / 2, 
            0
            ]
        self.left_mesh = Cylinder(vertical_length[0], inner_radius[0],
                                  position=left_mesh_position)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            -vertical_separation[0] / 2, 
            vertical_length[0] / 2, 
            0
            ]
        self.right_mesh = Cylinder(vertical_length[0], inner_radius[0],
                                  position=right_mesh_position)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        curve_mesh_position = [
            0, 
            vertical_length[0], 
            0
            ]
        curve_mesh_rotation = [-np.pi / 2, 0, 0]
        self.curve_mesh = Torus(vertical_separation[0] / 2, inner_radius[0], np.pi,
                                position=curve_mesh_position,
                                rotation=curve_mesh_rotation)
        vertices_list.append(self.curve_mesh.vertices)
        faces_list.append(self.curve_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.curve_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Flat_U_Handle(ConceptTemplate):
    def __init__(self, vertical_size, vertical_separation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.vertical_size = vertical_size
        self.vertical_separation = vertical_separation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        left_mesh_position = [
            vertical_separation[0] / 2, 
            vertical_size[1] / 2, 
            0
            ]
        self.left_mesh = Cuboid(vertical_size[1], vertical_size[0], vertical_size[2], 
                                position=left_mesh_position)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            -vertical_separation[0] / 2, 
            vertical_size[1] / 2, 
            0
            ]
        self.right_mesh = Cuboid(vertical_size[1], vertical_size[0], vertical_size[2], 
                                position=right_mesh_position)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        outer_radius = (vertical_separation[0] + vertical_size[0]) / 2
        inner_radius = (vertical_separation[0] - vertical_size[0]) / 2
        curve_mesh_position = [
            0, 
            vertical_size[1], 
            0
            ]
        curve_mesh_rotation = [-np.pi / 2, 0, 0]
        self.curve_mesh = Ring(vertical_size[2], outer_radius, inner_radius, np.pi,
                               position=curve_mesh_position,
                               rotation=curve_mesh_rotation)
        vertices_list.append(self.curve_mesh.vertices)
        faces_list.append(self.curve_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.curve_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Single_Cylinder(ConceptTemplate):
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

        self.mesh = Cylinder(size[1], size[0])
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Cylinder'