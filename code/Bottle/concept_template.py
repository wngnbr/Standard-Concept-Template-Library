import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh

class Multilevel_Body(ConceptTemplate):
    def __init__(self, num_levels, level_1_size, level_2_size, level_3_size, level_4_size, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.num_levels = num_levels
        self.level_1_size = level_1_size
        self.level_2_size = level_2_size
        self.level_3_size = level_3_size
        self.level_4_size = level_4_size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.bottom_mesh = Cylinder(level_1_size[1], level_1_size[0], level_1_size[2])
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        delta_height = level_1_size[1] / 2
        for i in range(num_levels[0] - 1):
            delta_height += locals()['level_'+ str(i+2) +'_size'][1] / 2
            mesh_position = [0, delta_height, 0]
            delta_height += locals()['level_'+ str(i+2) +'_size'][1] / 2
            self.mesh = Cylinder(locals()['level_'+ str(i+2) +'_size'][1], locals()['level_'+ str(i+2) +'_size'][0], locals()['level_'+ str(i+1) +'_size'][0], 
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

        self.semantic = 'Body'


class Cylindrical_Lid(ConceptTemplate):
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

        middle_radius = outer_size[1] * (1 - inner_size[2] / outer_size[2]) + outer_size[0] * inner_size[2] / outer_size[2]
        top_height = outer_size[2] - inner_size[2]

        bottom_mesh_position = [0, -(outer_size[2] - inner_size[2]) / 2, 0]
        self.bottom_mesh = Ring(inner_size[2], middle_radius, inner_size[0], 
                                outer_bottom_radius = outer_size[1],
                                inner_bottom_radius = inner_size[1],
                                position=bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [0, inner_size[2] / 2, 0]
        self.top_mesh = Cylinder(top_height, outer_size[0], middle_radius,
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

        self.semantic = 'Lid'