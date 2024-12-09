import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh

class Cylindrical_Barrel(ConceptTemplate):
    def __init__(self, size, thickness, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.thickness = thickness

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.mesh = Ring(size[2], size[0], size[0] - thickness[0],
                         outer_bottom_radius = size[1],
                         inner_bottom_radius = size[1] - thickness[0])
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


class Double_Layer_Barrel(ConceptTemplate):
    def __init__(self, main_size, bottom_size, thickness, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.main_size = main_size
        self.bottom_size = bottom_size
        self.thickness = thickness

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.mesh = Ring(main_size[2], main_size[0], main_size[0] - thickness[0],
                         outer_bottom_radius = main_size[1],
                         inner_bottom_radius = main_size[1] - thickness[0])
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        bottom_mesh_position = [
            0,
            -(main_size[2] + bottom_size[1]) / 2,
            0
        ]
        self.bottom_mesh = Ring(bottom_size[1], main_size[1], main_size[1] - thickness[0],
                                outer_bottom_radius = bottom_size[0],
                                inner_bottom_radius = bottom_size[0] - thickness[0],
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


class Single_Cap(ConceptTemplate):
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

        self.semantic = 'Cap'


class Trifold_Clip(ConceptTemplate):
    def __init__(self, clip_root_size, clip_vertical_size, clip_tip_size, clip_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.clip_root_size = clip_root_size
        self.clip_vertical_size = clip_vertical_size
        self.clip_tip_size = clip_tip_size
        self.clip_offset = clip_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0,
            -clip_root_size[1] / 2,
            clip_root_size[2] / 2
        ]
        self.mesh = Cuboid(clip_root_size[1], clip_root_size[0], clip_root_size[2],
                           position = mesh_position)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        middle_mesh_position = [
            0,
            -clip_vertical_size[1] / 2 - clip_offset[0],
            clip_root_size[2] + clip_vertical_size[2] / 2
        ]
        self.middle_mesh = Cuboid(clip_vertical_size[1], clip_vertical_size[0], clip_vertical_size[2],
                                  position = middle_mesh_position)
        vertices_list.append(self.middle_mesh.vertices)
        faces_list.append(self.middle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.middle_mesh.vertices)

        tip_mesh_position = [
            0,
            -clip_vertical_size[1] - clip_offset[0] + clip_tip_size[1] / 2 - clip_offset[1],
            clip_root_size[2] - clip_tip_size[2] / 2
        ]
        self.tip_mesh = Cuboid(clip_tip_size[1], clip_tip_size[0], clip_tip_size[2],
                               position = tip_mesh_position)
        vertices_list.append(self.tip_mesh.vertices)
        faces_list.append(self.tip_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.tip_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Clip'


class Curved_Clip(ConceptTemplate):
    def __init__(self, clip_curve_size, clip_curve_exist_angle, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        clip_curve_exist_angle = [x / 180 * np.pi for x in clip_curve_exist_angle]
        super().__init__(position, rotation)

        # Record Parameters
        self.clip_curve_size = clip_curve_size
        self.clip_curve_exist_angle = clip_curve_exist_angle

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0, 
            -clip_curve_size[0] * np.sin(clip_curve_exist_angle[0] / 2), 
            -clip_curve_size[0] * np.cos(clip_curve_exist_angle[0] / 2)
        ]
        mesh_rotation = [
            0,
            -np.pi / 2 + clip_curve_exist_angle[0] / 2,
            np.pi / 2
        ]
        self.mesh = Ring(clip_curve_size[2], clip_curve_size[0], clip_curve_size[1], clip_curve_exist_angle[0],
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

        self.semantic = 'Clip'


class Cylindrical_Button(ConceptTemplate):
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
            size[2] / 2,
            0
        ]
        self.mesh = Cylinder(size[2], size[0], size[1],
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

        self.semantic = 'Button'


class Bistratal_Button(ConceptTemplate):
    def __init__(self, bottom_size, top_size, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.bottom_size = bottom_size
        self.top_size = top_size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0,
            bottom_size[2] / 2,
            0
        ]
        self.mesh = Cylinder(bottom_size[2], bottom_size[0], bottom_size[1],
                             position = mesh_position)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        top_mesh_position = [
            0,
            bottom_size[2] + top_size[2] / 2,
            0
        ]
        self.top_mesh = Cylinder(top_size[2], top_size[0], top_size[1],
                                 position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Button'


class Cylindrical_Refill(ConceptTemplate):
    def __init__(self, size, tip_radius, tip_height, tip_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.tip_radius = tip_radius
        self.tip_height = tip_height
        self.tip_offset = tip_offset

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

        bottom_mesh_position = [
            0,
            -size[1] / 2 - tip_height[0] / 2, 
            0
        ]
        self.bottom_mesh = Cylinder(tip_height[0], tip_radius[0], 
                                    position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        tip_mesh_position = [
            0,
            -size[1] / 2 - tip_height[0], 
            0
        ]
        tip_mesh_rotation = [0, 0, np.pi]
        self.tip_mesh = Cone(tip_radius[0], tip_height[1], tip_offset,
                             position = tip_mesh_position,
                             rotation = tip_mesh_rotation)
        vertices_list.append(self.tip_mesh.vertices)
        faces_list.append(self.tip_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.tip_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Refill'