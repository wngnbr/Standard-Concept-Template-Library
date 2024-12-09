import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation, adjust_position_from_rotation, list_add
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


class Cylindrical_Cover(ConceptTemplate):
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

        self.semantic = 'Cover'


class Carved_Cylindrical_Cover(ConceptTemplate):
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
        top_mesh_position = [0, inner_size[2] / 2, 0]
        self.top_mesh = Cylinder(top_height, outer_size[0], middle_radius,
                                 position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [0, -(outer_size[2] - inner_size[2]) / 2, 0]
        self.bottom_mesh = Ring(inner_size[2], middle_radius, inner_size[0], 
                             outer_bottom_radius = outer_size[1],
                             inner_bottom_radius = inner_size[1],
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


class Semi_Spherical_Cover(ConceptTemplate):
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

        mesh_position = [
            0,
            -radius[0] * np.cos(exist_angle[0]),
            0
        ]
        self.mesh = Sphere(radius[0], 0, exist_angle[0],
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


class Cuboidal_Tophandle(ConceptTemplate):
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
            0
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

        self.semantic = 'Handle'


class Trifold_Tophandle(ConceptTemplate):
    def __init__(self, mounting_size, mounting_seperation, grip_size, mounting_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        mounting_rotation = [x / 180 * np.pi for x in mounting_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.mounting_size = mounting_size
        self.mounting_seperation = mounting_seperation
        self.grip_size = grip_size
        self.mounting_rotation = mounting_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        left_mesh_position = [
            mounting_seperation[0] / 2 + mounting_size[1] * np.sin(mounting_rotation[0]) / 2, 
            mounting_size[1] * np.cos(mounting_rotation[0]) / 2, 
            0
        ]
        left_mesh_rotation = [0, 0, -mounting_rotation[0]]
        self.left_mesh = Cuboid(mounting_size[1], mounting_size[0], mounting_size[2], 
                                position = left_mesh_position,
                                rotation = left_mesh_rotation)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            -mounting_seperation[0] / 2 - mounting_size[1] * np.sin(mounting_rotation[0]) / 2, 
            mounting_size[1] * np.cos(mounting_rotation[0]) / 2, 
            0
        ]
        right_mesh_rotation = [0, 0, mounting_rotation[0]]
        self.right_mesh = Cuboid(mounting_size[1], mounting_size[0], mounting_size[2], 
                                 position = right_mesh_position,
                                 rotation = right_mesh_rotation)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        top_mesh_position = [
            0, 
            mounting_size[1] * np.cos(mounting_rotation[0]) + grip_size[1] / 2, 
            0
        ]
        self.top_mesh = Cuboid(grip_size[1], grip_size[0], grip_size[2], 
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

        self.semantic = 'Handle'


class Semi_Ring_Tophandle(ConceptTemplate):
    def __init__(self, curve_size, curve_exist_angle, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        curve_exist_angle = [x / 180 * np.pi for x in curve_exist_angle]
        super().__init__(position, rotation)

        # Record Parameters
        self.curve_size = curve_size
        self.curve_exist_angle = curve_exist_angle

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0,
            -curve_size[0] * np.cos(curve_exist_angle[0] / 2),
            0
        ]
        mesh_rotation = [
            -np.pi / 2,
            -np.pi / 2 + curve_exist_angle[0] / 2,
            0
        ]
        self.mesh = Ring(curve_size[2], curve_size[0], curve_size[1], curve_exist_angle[0],
                         position = mesh_position,
                         rotation = mesh_rotation,
                         rotation_order = "YXZ")
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


class Multilevel_Tophandle(ConceptTemplate):
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

        delta_height = 0
        for i in range(num_levels[0]):
            delta_height += locals()['level_'+ str(i+1) +'_size'][2] / 2
            mesh_position = [0, delta_height, 0]
            delta_height += locals()['level_'+ str(i+1) +'_size'][2] / 2
            self.mesh = Cylinder(locals()['level_'+ str(i+1) +'_size'][2], locals()['level_'+ str(i+1) +'_size'][0], locals()['level_'+ str(i+1) +'_size'][1], 
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

        self.semantic = 'Handle'


class Trifold_Sidehandle(ConceptTemplate):
    def __init__(self, mounting_size, mounting_seperation, grip_size, mounting_rotation, handle_seperation, whole_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        mounting_rotation = [x / 180 * np.pi for x in mounting_rotation]
        whole_rotation = [x / 180 * np.pi for x in whole_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.mounting_size = mounting_size
        self.mounting_seperation = mounting_seperation
        self.grip_size = grip_size
        self.mounting_rotation = mounting_rotation
        self.handle_seperation = handle_seperation
        self.whole_rotation = whole_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # level_1
        level_1_position = [handle_seperation[0] / 2, 0, 0]
        level_1_rotation = [0, 0, whole_rotation[0]]

        left_1_mesh_position = [
            mounting_size[0] * np.cos(mounting_rotation[0]) / 2, 
            0,
            mounting_seperation[0] / 2 + mounting_size[0] * np.sin(mounting_rotation[0]) / 2
        ]
        left_1_mesh_rotation = [0, -mounting_rotation[0], 0]

        left_1_mesh_position = adjust_position_from_rotation(left_1_mesh_position, level_1_rotation)
        left_1_mesh_position = list_add(left_1_mesh_position, level_1_position)
        left_1_mesh_rotation = list_add(left_1_mesh_rotation, level_1_rotation)

        self.left_1_mesh = Cuboid(mounting_size[1], mounting_size[0], mounting_size[2], 
                                  position = left_1_mesh_position,
                                  rotation = left_1_mesh_rotation)
        vertices_list.append(self.left_1_mesh.vertices)
        faces_list.append(self.left_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_1_mesh.vertices)

        right_1_mesh_position = [
            mounting_size[0] * np.cos(mounting_rotation[0]) / 2, 
            0, 
            -mounting_seperation[0] / 2 - mounting_size[0] * np.sin(mounting_rotation[0]) / 2, 
        ]
        right_1_mesh_rotation = [0, mounting_rotation[0], 0]

        right_1_mesh_position = adjust_position_from_rotation(right_1_mesh_position, level_1_rotation)
        right_1_mesh_position = list_add(right_1_mesh_position, level_1_position)
        right_1_mesh_rotation = list_add(right_1_mesh_rotation, level_1_rotation)

        self.right_1_mesh = Cuboid(mounting_size[1], mounting_size[0], mounting_size[2], 
                                 position = right_1_mesh_position,
                                 rotation = right_1_mesh_rotation)
        vertices_list.append(self.right_1_mesh.vertices)
        faces_list.append(self.right_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_1_mesh.vertices)

        top_1_mesh_position = [
            mounting_size[0] * np.cos(mounting_rotation[0]) + grip_size[0] / 2, 
            0, 
            0
        ]
        top_1_mesh_rotation = [0, 0, 0]

        top_1_mesh_position = adjust_position_from_rotation(top_1_mesh_position, level_1_rotation)
        top_1_mesh_position = list_add(top_1_mesh_position, level_1_position)
        top_1_mesh_rotation = list_add(top_1_mesh_rotation, level_1_rotation)

        self.top_1_mesh = Cuboid(grip_size[1], grip_size[0], grip_size[2], 
                                 position = top_1_mesh_position,
                                 rotation = top_1_mesh_rotation)
        vertices_list.append(self.top_1_mesh.vertices)
        faces_list.append(self.top_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_1_mesh.vertices)

        # level_2
        level_2_position = [-handle_seperation[0] / 2, 0, 0]
        level_2_rotation = [0, 0, -whole_rotation[0]]

        left_2_mesh_position = [
            -mounting_size[0] * np.cos(mounting_rotation[0]) / 2, 
            0,
            mounting_seperation[0] / 2 + mounting_size[0] * np.sin(mounting_rotation[0]) / 2
        ]
        left_2_mesh_rotation = [0, mounting_rotation[0], 0]

        left_2_mesh_position = adjust_position_from_rotation(left_2_mesh_position, level_2_rotation)
        left_2_mesh_position = list_add(left_2_mesh_position, level_2_position)
        left_2_mesh_rotation = list_add(left_2_mesh_rotation, level_2_rotation)

        self.left_2_mesh = Cuboid(mounting_size[1], mounting_size[0], mounting_size[2], 
                                  position = left_2_mesh_position,
                                  rotation = left_2_mesh_rotation)
        vertices_list.append(self.left_2_mesh.vertices)
        faces_list.append(self.left_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_2_mesh.vertices)

        right_2_mesh_position = [
            -mounting_size[0] * np.cos(mounting_rotation[0]) / 2, 
            0, 
            -mounting_seperation[0] / 2 - mounting_size[0] * np.sin(mounting_rotation[0]) / 2, 
        ]
        right_2_mesh_rotation = [0, -mounting_rotation[0], 0]

        right_2_mesh_position = adjust_position_from_rotation(right_2_mesh_position, level_2_rotation)
        right_2_mesh_position = list_add(right_2_mesh_position, level_2_position)
        right_2_mesh_rotation = list_add(right_2_mesh_rotation, level_2_rotation)

        self.right_2_mesh = Cuboid(mounting_size[1], mounting_size[0], mounting_size[2], 
                                 position = right_2_mesh_position,
                                 rotation = right_2_mesh_rotation)
        vertices_list.append(self.right_2_mesh.vertices)
        faces_list.append(self.right_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_2_mesh.vertices)

        top_2_mesh_position = [
            -mounting_size[0] * np.cos(mounting_rotation[0]) - grip_size[0] / 2, 
            0, 
            0
        ]
        top_2_mesh_rotation = [0, 0, 0]

        top_2_mesh_position = adjust_position_from_rotation(top_2_mesh_position, level_2_rotation)
        top_2_mesh_position = list_add(top_2_mesh_position, level_2_position)
        top_2_mesh_rotation = list_add(top_2_mesh_rotation, level_2_rotation)

        self.top_2_mesh = Cuboid(grip_size[1], grip_size[0], grip_size[2], 
                                 position = top_2_mesh_position,
                                 rotation = top_2_mesh_rotation)
        vertices_list.append(self.top_2_mesh.vertices)
        faces_list.append(self.top_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_2_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class L_Shaped_Sidehandle(ConceptTemplate):
    def __init__(self, bottom_size, top_size, handle_seperation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.bottom_size = bottom_size
        self.top_size = top_size
        self.handle_seperation = handle_seperation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # level_1
        bottom_1_mesh_position = [
            handle_seperation[0] / 2 + bottom_size[0] / 2, 
            0,
            0
        ]
        self.bottom_1_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2], 
                                    position = bottom_1_mesh_position)
        vertices_list.append(self.bottom_1_mesh.vertices)
        faces_list.append(self.bottom_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_1_mesh.vertices)

        top_1_mesh_position = [
            handle_seperation[0] / 2 + top_size[0] / 2, 
            (bottom_size[1] + top_size[1]) / 2,
            0
        ]

        self.top_1_mesh = Cuboid(top_size[1], top_size[0], top_size[2], 
                                 position = top_1_mesh_position)
        vertices_list.append(self.top_1_mesh.vertices)
        faces_list.append(self.top_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_1_mesh.vertices)

        # level_2
        bottom_2_mesh_position = [
            -handle_seperation[0] / 2 - bottom_size[0] / 2, 
            0,
            0
        ]
        self.bottom_2_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2], 
                                    position = bottom_2_mesh_position)
        vertices_list.append(self.bottom_2_mesh.vertices)
        faces_list.append(self.bottom_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_2_mesh.vertices)

        top_2_mesh_position = [
            -handle_seperation[0] / 2 - top_size[0] / 2, 
            (bottom_size[1] + top_size[1]) / 2,
            0
        ]

        self.top_2_mesh = Cuboid(top_size[1], top_size[0], top_size[2], 
                                 position = top_2_mesh_position)
        vertices_list.append(self.top_2_mesh.vertices)
        faces_list.append(self.top_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_2_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Cuboidal_Sidehandle(ConceptTemplate):
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
            0,
            -size[2] / 2
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

        self.semantic = 'Handle'