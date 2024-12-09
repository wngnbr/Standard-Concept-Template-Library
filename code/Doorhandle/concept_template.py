import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh


class Regular_handle(ConceptTemplate):
    def __init__(self, fixed_part_size, vertical_movable_size, horizontal_movable_size, interpiece_offset_1, interpiece_offset_2, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.fixed_part_size = fixed_part_size
        self.vertical_movable_size = vertical_movable_size
        self.horizontal_movable_size = horizontal_movable_size
        self.interpiece_offset_1 = interpiece_offset_1
        self.interpiece_offset_2 = interpiece_offset_2

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        base_mesh_position = [0, 0, fixed_part_size[2] / 2]
        self.base_mesh = Cuboid(fixed_part_size[1], fixed_part_size[0], fixed_part_size[2],
                                position=base_mesh_position)
        vertices_list.append(self.base_mesh.vertices)
        faces_list.append(self.base_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.base_mesh.vertices)

        middle_mesh_position = [
            interpiece_offset_1[0],
            interpiece_offset_1[1],
            fixed_part_size[2] + vertical_movable_size[2] / 2,
        ]
        self.middle_mesh = Cuboid(vertical_movable_size[1], vertical_movable_size[0], vertical_movable_size[2],
                                  position=middle_mesh_position)
        vertices_list.append(self.middle_mesh.vertices)
        faces_list.append(self.middle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.middle_mesh.vertices)

        main_mesh_position = [
            interpiece_offset_1[0] + interpiece_offset_2[0],
            interpiece_offset_1[1] + interpiece_offset_2[1],
            fixed_part_size[2] + vertical_movable_size[2] + horizontal_movable_size[2] / 2,
        ]
        self.main_mesh = Cuboid(horizontal_movable_size[1], horizontal_movable_size[0], horizontal_movable_size[2],
                                position=main_mesh_position)
        vertices_list.append(self.main_mesh.vertices)
        faces_list.append(self.main_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.main_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Knob_handle(ConceptTemplate):
    def __init__(self, fixed_part_size, sub_size, main_size, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.fixed_part_size = fixed_part_size
        self.sub_size = sub_size
        self.main_size = main_size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        base_mesh_position = [0, 0, fixed_part_size[1] / 2]
        base_mesh_rotation = [np.pi / 2, 0, 0]
        self.base_mesh = Cylinder(fixed_part_size[1], fixed_part_size[0], fixed_part_size[0],
                                  position=base_mesh_position,
                                  rotation=base_mesh_rotation)
        vertices_list.append(self.base_mesh.vertices)
        faces_list.append(self.base_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.base_mesh.vertices)

        middle_mesh_position = [0, 0, fixed_part_size[1] + sub_size[1] / 2]
        middle_mesh_rotation = [np.pi / 2, 0, 0]
        self.middle_mesh = Cylinder(sub_size[1], sub_size[0], sub_size[0],
                                    position=middle_mesh_position,
                                    rotation=middle_mesh_rotation)
        vertices_list.append(self.middle_mesh.vertices)
        faces_list.append(self.middle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.middle_mesh.vertices)

        main_mesh_position = [0, 0, fixed_part_size[1] + sub_size[1] + main_size[1] / 2]
        main_mesh_rotation = [np.pi / 2, 0, 0]
        self.main_mesh = Cylinder(main_size[1], main_size[0], main_size[0],
                                  position=main_mesh_position,
                                  rotation=main_mesh_rotation)
        vertices_list.append(self.main_mesh.vertices)
        faces_list.append(self.main_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.main_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class TShaped_handle(ConceptTemplate):
    def __init__(self, sub_size, main_size, interpiece_offset, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.sub_size = sub_size
        self.main_size = main_size
        self.interpiece_offset = interpiece_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        sub_mesh_position = [0, 0, sub_size[1] / 2]
        sub_mesh_rotation = [np.pi / 2, 0, 0]
        self.sub_mesh = Cylinder(sub_size[1], sub_size[0] / 2, sub_size[0] / 2,
                                 rotation=sub_mesh_rotation,
                                 position=sub_mesh_position)
        vertices_list.append(self.sub_mesh.vertices)
        faces_list.append(self.sub_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.sub_mesh.vertices)

        main_mesh_position = [
            interpiece_offset[0],
            interpiece_offset[1],
            sub_size[1] + main_size[2] / 2,
        ]
        self.main_mesh = Cuboid(main_size[1], main_size[0], main_size[2],
                                position=main_mesh_position)
        vertices_list.append(self.main_mesh.vertices)
        faces_list.append(self.main_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.main_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'