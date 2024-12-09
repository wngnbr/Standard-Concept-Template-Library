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


class Multilevel_Body(ConceptTemplate):
    def __init__(self, num_levels, level_1_bottom_radius, level_1_top_radius, level_1_height, level_2_top_radius=0, level_2_height=0, level_3_top_radius=0, level_3_height=0, level_4_top_radius=0, level_4_height=0, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.num_levels = num_levels
        self.level_1_bottom_radius = level_1_bottom_radius
        self.level_1_top_radius = level_1_top_radius
        self.level_1_height = level_1_height
        self.level_2_top_radius = level_2_top_radius
        self.level_2_height = level_2_height
        self.level_3_top_radius = level_3_top_radius
        self.level_3_height = level_3_height
        self.level_4_top_radius = level_4_top_radius
        self.level_4_height = level_4_height

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_1_top_radius = level_1_top_radius[0] * (1 - level_1_height[1] / level_1_height[0]) + level_1_bottom_radius[0] * level_1_height[1] / level_1_height[0]
        mesh_1_height = level_1_height[0] - level_1_height[1]
        bottom_mesh_position = [0, -level_1_height[1] / 2, 0]
        self.bottom_mesh = Cylinder(mesh_1_height, mesh_1_top_radius, level_1_bottom_radius[0], 
                                    position=bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [0, (level_1_height[0] - level_1_height[1]) / 2, 0]
        self.top_mesh = Ring(level_1_height[1], level_1_top_radius[0], level_1_top_radius[1], 
                             outer_bottom_radius = mesh_1_top_radius,
                             inner_bottom_radius = level_1_bottom_radius[1],
                             position=top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        delta_height = level_1_height[0] / 2
        for i in range(num_levels[0] - 1):
            delta_height += locals()['level_'+ str(i+2) +'_height'][0] / 2
            top_mesh_position = [0, delta_height, 0]
            delta_height += locals()['level_'+ str(i+2) +'_height'][0] / 2
            self.top_mesh = Ring(locals()['level_'+ str(i+2) +'_height'][0], locals()['level_'+ str(i+2) +'_top_radius'][0], locals()['level_'+ str(i+2) +'_top_radius'][1], 
                                 outer_bottom_radius = locals()['level_'+ str(i+1) +'_top_radius'][0],
                                 inner_bottom_radius = locals()['level_'+ str(i+1) +'_top_radius'][1],
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
    def __init__(self, horizontal_thickness, horizontal_length, vertical_thickness, horizontal_rotation, horizontal_separation, mounting_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        horizontal_rotation = [x / 180 * np.pi for x in horizontal_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.horizontal_thickness = horizontal_thickness
        self.horizontal_length = horizontal_length
        self.vertical_thickness = vertical_thickness
        self.horizontal_rotation = horizontal_rotation
        self.horizontal_separation = horizontal_separation
        self.mounting_offset = mounting_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        top_mesh_position = [
            0, 
            horizontal_separation[0] / 2 - horizontal_length[0] * np.sin(horizontal_rotation[0]) / 2, 
            mounting_offset[0] + horizontal_length[0] * np.cos(horizontal_rotation[0]) / 2
            ]
        top_mesh_rotation = [horizontal_rotation[0], 0, 0]
        self.top_mesh = Cuboid(horizontal_thickness[1], horizontal_thickness[0], horizontal_length[0], 
                               position=top_mesh_position,
                               rotation=top_mesh_rotation)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0, 
            -horizontal_separation[0] / 2 - horizontal_length[1] * np.sin(horizontal_rotation[1]) / 2, 
            horizontal_length[1] * np.cos(horizontal_rotation[1]) / 2
            ]
        bottom_mesh_rotation = [horizontal_rotation[1], 0, 0]
        self.bottom_mesh = Cuboid(horizontal_thickness[1], horizontal_thickness[0], horizontal_length[1], 
                                  position=bottom_mesh_position,
                                  rotation=bottom_mesh_rotation)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        delta_y = horizontal_separation[0] - horizontal_length[0] * np.sin(horizontal_rotation[0]) + horizontal_length[1] * np.sin(horizontal_rotation[1])
        delta_z = mounting_offset[0] - horizontal_length[1] * np.cos(horizontal_rotation[1]) + horizontal_length[0] * np.cos(horizontal_rotation[0])
        vertical_length = np.sqrt(delta_y * delta_y + delta_z * delta_z) + horizontal_thickness[1]
        vertical_rotation = np.arctan(delta_z / delta_y)
        vertical_y_offset = (-horizontal_length[0] * np.sin(horizontal_rotation[0]) - horizontal_length[1] * np.sin(horizontal_rotation[1])) / 2
        vertical_z_offset = (horizontal_length[1] * np.cos(horizontal_rotation[1]) + mounting_offset[0] + horizontal_length[0] * np.cos(horizontal_rotation[0])) / 2
        vertical_mesh_position = [
            0, 
            vertical_y_offset, 
            vertical_z_offset + vertical_thickness[1] / 2
            ]
        vertical_mesh_rotation = [vertical_rotation, 0, 0]
        self.vertical_mesh = Cuboid(vertical_length, vertical_thickness[0], vertical_thickness[1], 
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
    def __init__(self, radius, central_angle, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        central_angle = [x / 180 * np.pi for x in central_angle]
        super().__init__(position, rotation)

        # Record Parameters
        self.radius = radius
        self.central_angle = central_angle

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_rotation = [0, 0, np.pi / 2]
        self.mesh = Torus(radius[0], radius[1], central_angle[0],
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