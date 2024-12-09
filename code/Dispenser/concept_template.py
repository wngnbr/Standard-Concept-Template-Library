import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation, adjust_position_from_rotation, list_add
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
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order='YXZ')

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Cuboidal_Body(ConceptTemplate):
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

        self.mesh = Cuboid(size[1], size[0], size[2])
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order='YXZ')

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Press_Nozzle(ConceptTemplate):
    def __init__(self, num_levels, level_1_size, level_2_size, level_3_size, level_4_size, level_5_size, num_nozzles, nozzle_size, nozzle_length, nozzle_offset, nozzle_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        nozzle_rotation = [x / 180 * np.pi for x in nozzle_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.num_levels = num_levels
        self.level_1_size = level_1_size
        self.level_2_size = level_2_size
        self.level_3_size = level_3_size
        self.level_4_size = level_4_size
        self.level_5_size = level_5_size
        self.num_nozzles = num_nozzles
        self.nozzle_size = nozzle_size
        self.nozzle_length = nozzle_length
        self.nozzle_offset = nozzle_offset
        self.nozzle_rotation = nozzle_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        delta_height = 0
        for i in range(num_levels[0]):
            delta_height += locals()['level_'+ str(i+1) +'_size'][1] / 2
            mesh_position = [0, delta_height, 0]
            delta_height += locals()['level_'+ str(i+1) +'_size'][1] / 2
            self.mesh = Cylinder(locals()['level_'+ str(i+1) +'_size'][1], locals()['level_'+ str(i+1) +'_size'][0], 
                                 position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        total_offset_z = nozzle_length[0] / 2 * np.cos(nozzle_rotation[0]) + locals()['level_%d_size'%(num_levels[0])][0]
        nozzle_mesh_position = [
            0, 
            delta_height + nozzle_offset[0] - nozzle_length[0] / 2 * np.sin(nozzle_rotation[0]), 
            total_offset_z
        ]
        nozzle_mesh_rotation = [nozzle_rotation[0], 0, 0]
        self.nozzle_mesh = Cuboid(nozzle_size[1], nozzle_size[0], nozzle_length[0],
                                  position = nozzle_mesh_position,
                                  rotation = nozzle_mesh_rotation)
        vertices_list.append(self.nozzle_mesh.vertices)
        faces_list.append(self.nozzle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.nozzle_mesh.vertices)

        if num_nozzles[0] == 2:
            total2_offset_z = nozzle_length[0] * np.cos(nozzle_rotation[0]) + locals()['level_%d_size'%(num_levels[0])][0] + nozzle_length[1] * np.cos(nozzle_rotation[1]) / 2
            nozzle_mesh_position = [
                0, 
                delta_height + nozzle_offset[0] - nozzle_length[0] * np.sin(nozzle_rotation[0]) - nozzle_length[1] / 2 * np.sin(nozzle_rotation[1]), 
                total2_offset_z
            ]
            nozzle_mesh_rotation = [nozzle_rotation[1], 0, 0]
            self.nozzle_mesh = Cuboid(nozzle_size[1], nozzle_size[0], nozzle_length[1],
                                      position = nozzle_mesh_position,
                                      rotation = nozzle_mesh_rotation)
            vertices_list.append(self.nozzle_mesh.vertices)
            faces_list.append(self.nozzle_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.nozzle_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order='YXZ')

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Nozzle'


class Spray_Nozzle(ConceptTemplate):
    def __init__(self, bottom_size, middle_size, top_size, top_offset, top_rotation, nozzle_size, handle_size, handle_offset, handle_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        top_rotation = [x / 180 * np.pi for x in top_rotation]
        handle_rotation = [x / 180 * np.pi for x in handle_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.bottom_size = bottom_size
        self.middle_size = middle_size
        self.top_size = top_size
        self.top_offset = top_offset
        self.top_rotation = top_rotation
        self.nozzle_size = nozzle_size
        self.handle_size = handle_size
        self.handle_offset = handle_offset
        self.handle_rotation = handle_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # bottom
        delta_height = bottom_size[1] / 2
        bottom_mesh_position = [
            0, 
            delta_height,
            0
        ]
        self.bottom_mesh = Cylinder(bottom_size[1], bottom_size[0], 
                                    position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        # middle
        delta_height += bottom_size[1] / 2 + middle_size[1] / 2
        middle_mesh_position = [
            0, 
            delta_height,
            0
        ]
        self.middle_mesh = Cuboid(middle_size[1], middle_size[0], middle_size[2],
                                  position = middle_mesh_position)
        vertices_list.append(self.middle_mesh.vertices)
        faces_list.append(self.middle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.middle_mesh.vertices)

        # top
        top_offset_y = -top_offset[0] * np.tan(top_rotation[0])
        delta_height += middle_size[1] / 2 + top_size[1] / 2
        top_mesh_position = [
            0, 
            delta_height + top_offset_y + top_offset[0] * np.sin(top_rotation[0]),
            top_offset[0] * np.cos(top_rotation[0])
        ]
        top_mesh_rotation = [top_rotation[0], 0, 0]
        self.top_mesh = Cuboid(top_size[1], top_size[0], top_size[2],
                               position = top_mesh_position,
                               rotation = top_mesh_rotation)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        # nozzle
        nozzle_total_offset_z = top_offset[0] + (top_size[2] + nozzle_size[1]) / 2
        nozzle_mesh_position = [
            0, 
            delta_height + top_offset_y - nozzle_total_offset_z * np.sin(top_rotation[0]),
            nozzle_total_offset_z * np.cos(top_rotation[0])
        ]
        nozzle_mesh_rotation = [top_rotation[0] + np.pi / 2, 0, 0]
        self.nozzle_mesh = Cylinder(nozzle_size[1], nozzle_size[0], 
                                    position = nozzle_mesh_position,
                                    rotation = nozzle_mesh_rotation)
        vertices_list.append(self.nozzle_mesh.vertices)
        faces_list.append(self.nozzle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.nozzle_mesh.vertices)

        # handle
        handle_total_offset_z = top_offset[0] + top_size[2] / 2 + handle_offset[0] - handle_size[2] / 2 - handle_size[1] / 2 * np.sin(handle_rotation[0])
        handle_total_offset_y = -handle_size[1] / 2 * np.cos(handle_rotation[0]) - top_size[1] / 2
        handle_mesh_position_1 = [
            0, 
            handle_total_offset_y,
            handle_total_offset_z
        ]
        handle_mesh_rotation = [top_rotation[0], 0, 0]
        handle_mesh_position_1 = adjust_position_from_rotation(handle_mesh_position_1, handle_mesh_rotation)

        handle_mesh_position_2 = [
            0, 
            delta_height + top_offset_y,
            0
        ]
        handle_mesh_position = list_add(handle_mesh_position_1, handle_mesh_position_2)
        handle_mesh_rotation[0] += handle_rotation[0]

        self.handle_mesh = Cuboid(handle_size[1], handle_size[0], handle_size[2],
                                  position = handle_mesh_position,
                                  rotation = handle_mesh_rotation)
        vertices_list.append(self.handle_mesh.vertices)
        faces_list.append(self.handle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.handle_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order='YXZ')

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Nozzle'