import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh

class Semi_Spherical_Body(ConceptTemplate):
    def __init__(self, horizontal_axis, vertical_axis, exist_angle, bottom_size, thickness, x_z_ratio, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        exist_angle = [x / 180 * np.pi for x in exist_angle]
        super().__init__(position, rotation)

        # Record Parameters
        self.horizontal_axis = horizontal_axis
        self.vertical_axis = vertical_axis
        self.exist_angle = exist_angle
        self.bottom_size = bottom_size
        self.thickness = thickness
        self.x_z_ratio = x_z_ratio

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.top_mesh = Sphere(horizontal_axis[0] * x_z_ratio[0], np.pi / 2 - exist_angle[0], np.pi / 2,
                               radius_y = vertical_axis[0], 
                               radius_z = horizontal_axis[0])
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        middle_mesh_rotation = [np.pi, 0, 0]
        self.middle_mesh = Sphere(horizontal_axis[0] * x_z_ratio[0], np.pi / 2 - exist_angle[1], np.pi / 2,
                                  radius_y = vertical_axis[1], 
                                  radius_z = horizontal_axis[0],
                                  rotation = middle_mesh_rotation)
        vertices_list.append(self.middle_mesh.vertices)
        faces_list.append(self.middle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.middle_mesh.vertices)

        bottom_radius = horizontal_axis[0] * np.cos(exist_angle[1])
        bottom_offset = horizontal_axis[0] * np.sin(exist_angle[1]) * vertical_axis[1] / horizontal_axis[0]
        bottom_mesh_position = [
            0, 
            -bottom_offset - bottom_size[1] / 2, 
            0
        ]
        self.bottom_mesh = Cylinder(bottom_size[1], bottom_radius * x_z_ratio[0], bottom_size[0] * x_z_ratio[0],
                                    top_radius_z = bottom_radius,
                                    bottom_radius_z = bottom_size[0],
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


class Spherical_Cylindrical_Body(ConceptTemplate):
    def __init__(self, horizontal_axis, vertical_axis, exist_angle, bottom_size, thickness, x_z_ratio, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        exist_angle = [x / 180 * np.pi for x in exist_angle]
        super().__init__(position, rotation)

        # Record Parameters
        self.horizontal_axis = horizontal_axis
        self.vertical_axis = vertical_axis
        self.exist_angle = exist_angle
        self.bottom_size = bottom_size
        self.thickness = thickness
        self.x_z_ratio = x_z_ratio

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.top_mesh = Sphere(horizontal_axis[0] * x_z_ratio[0], np.pi / 2 - exist_angle[0], np.pi / 2,
                               radius_y = vertical_axis[0], 
                               radius_z = horizontal_axis[0])
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0, 
            -bottom_size[1] / 2, 
            0
        ]
        self.bottom_mesh = Cylinder(bottom_size[1], horizontal_axis[0] * x_z_ratio[0], bottom_size[0] * x_z_ratio[0],
                                    top_radius_z = horizontal_axis[0],
                                    bottom_radius_z = bottom_size[0],
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


class Multilevel_Body(ConceptTemplate):
    def __init__(self, num_levels, level_1_bottom_radius, level_1_top_radius, level_1_height, level_2_top_radius=0, level_2_height=0, level_3_top_radius=0, level_3_height=0, level_4_top_radius=0, level_4_height=0, level_5_top_radius=0, level_5_height=0, x_z_ratio=1, position = [0, 0, 0], rotation = [0, 0, 0]):

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
        self.level_5_top_radius = level_5_top_radius
        self.level_5_height = level_5_height
        self.x_z_ratio = x_z_ratio

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_1_top_radius = level_1_top_radius[0] * (1 - level_1_height[1] / level_1_height[0]) + level_1_bottom_radius[0] * level_1_height[1] / level_1_height[0]
        mesh_1_height = level_1_height[0] - level_1_height[1]
        bottom_mesh_position = [0, -level_1_height[1] / 2, 0]
        self.bottom_mesh = Cylinder(mesh_1_height, mesh_1_top_radius * x_z_ratio[0], level_1_bottom_radius[0] * x_z_ratio[0], 
                                    top_radius_z = mesh_1_top_radius,
                                    bottom_radius_z = level_1_bottom_radius[0],
                                    position=bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [0, (level_1_height[0] - level_1_height[1]) / 2, 0]
        self.top_mesh = Ring(level_1_height[1], level_1_top_radius[0], level_1_top_radius[1], 
                             outer_bottom_radius = mesh_1_top_radius,
                             inner_bottom_radius = level_1_bottom_radius[1],
                             x_z_ratio = x_z_ratio[0],
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
                                 x_z_ratio = x_z_ratio[0],
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


class Standard_Cover(ConceptTemplate):
    def __init__(self, outer_size, inner_size, num_knobs, knob_1_size, knob_2_size, knob_3_size, knob_4_size, knob_5_size, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.outer_size = outer_size
        self.inner_size = inner_size
        self.num_knobs = num_knobs
        self.knob_1_size = knob_1_size
        self.knob_2_size = knob_2_size
        self.knob_3_size = knob_3_size
        self.knob_4_size = knob_4_size
        self.knob_5_size = knob_5_size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_1_bottom_radius = outer_size[1] * (1 - inner_size[2] / outer_size[2]) + outer_size[0] * inner_size[2] / outer_size[2]
        mesh_1_height = outer_size[2] - inner_size[2]
        top_mesh_position = [
            0, 
            inner_size[2] / 2 + outer_size[2] / 2, 
            0
        ]
        self.top_mesh = Cylinder(mesh_1_height, outer_size[0], mesh_1_bottom_radius, 
                                 position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0, 
            -(outer_size[2] - inner_size[2]) / 2 + outer_size[2] / 2, 
            0
        ]
        self.bottom_mesh = Ring(inner_size[2], mesh_1_bottom_radius, inner_size[0], 
                             outer_bottom_radius = outer_size[1],
                             inner_bottom_radius = inner_size[1],
                             position=bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        delta_height = outer_size[2]
        for i in range(num_knobs[0]):
            delta_height += locals()['knob_%d_size'%(i+1)][2] / 2
            knob_mesh_position = [0, delta_height, 0]
            self.knob_mesh = Cylinder(locals()['knob_%d_size'%(i+1)][2], locals()['knob_%d_size'%(i+1)][0], locals()['knob_%d_size'%(i+1)][1], 
                                      position=knob_mesh_position)
            delta_height += locals()['knob_%d_size'%(i+1)][2] / 2
            vertices_list.append(self.knob_mesh.vertices)
            faces_list.append(self.knob_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.knob_mesh.vertices)


        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Cover'


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
            -mounting_offset[0] - horizontal_length[0] * np.cos(horizontal_rotation[0]) / 2
            ]
        top_mesh_rotation = [-horizontal_rotation[0], 0, 0]
        self.top_mesh = Cuboid(horizontal_thickness[1], horizontal_thickness[0], horizontal_length[0], 
                               position=top_mesh_position,
                               rotation=top_mesh_rotation)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0, 
            -horizontal_separation[0] / 2 - horizontal_length[1] * np.sin(horizontal_rotation[1]) / 2, 
            -horizontal_length[1] * np.cos(horizontal_rotation[1]) / 2
            ]
        bottom_mesh_rotation = [-horizontal_rotation[1], 0, 0]
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
            -vertical_z_offset - vertical_thickness[1] / 2
            ]
        vertical_mesh_rotation = [-vertical_rotation, 0, 0]
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

        mesh_rotation = [np.pi, 0, -np.pi / 2]
        self.mesh = Torus(radius[0], radius[1], exist_angle[0],
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


class Ring_Handle(ConceptTemplate):
    def __init__(self, size, exist_angle, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        exist_angle = [x / 180 * np.pi for x in exist_angle]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.exist_angle = exist_angle

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_rotation = [np.pi, 0, -np.pi / 2]
        self.mesh = Ring(size[2], size[0], size[1], exist_angle[0],
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


class Cylindrical_Handle(ConceptTemplate):
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

        mesh_position = [0, 0, -size[2] / 2]
        mesh_rotation = [-np.pi / 2, 0, 0]
        self.mesh = Cylinder(size[2], size[0], size[1],
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

        self.semantic = 'Handle'


class Round_U_Handle(ConceptTemplate):
    def __init__(self, mounting_radius, vertical_separation, vertical_length, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.mounting_radius = mounting_radius
        self.vertical_separation = vertical_separation
        self.vertical_length = vertical_length

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        left_mesh_position = [
            0, 
            vertical_length[0] / 2, 
            vertical_separation[0] / 2
            ]
        self.left_mesh = Cylinder(vertical_length[0], mounting_radius[0],
                                  position=left_mesh_position)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            0, 
            vertical_length[0] / 2, 
            -vertical_separation[0] / 2
            ]
        self.right_mesh = Cylinder(vertical_length[0], mounting_radius[0],
                                  position=right_mesh_position)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        curve_mesh_position = [
            0, 
            vertical_length[0], 
            0
        ]
        curve_mesh_rotation = [0, np.pi / 2, 0]
        self.curve_mesh = Torus(vertical_separation[0] / 2, mounting_radius[0], np.pi,
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
            0, 
            vertical_size[1] / 2, 
            vertical_separation[0] / 2
            ]
        self.left_mesh = Cuboid(vertical_size[1], vertical_size[0], vertical_size[2], 
                                position=left_mesh_position)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            0, 
            vertical_size[1] / 2, 
            -vertical_separation[0] / 2
            ]
        self.right_mesh = Cuboid(vertical_size[1], vertical_size[0], vertical_size[2], 
                                position=right_mesh_position)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        outer_radius = (vertical_separation[0] + vertical_size[0]) / 2
        mounting_radius = (vertical_separation[0] - vertical_size[0]) / 2
        curve_mesh_position = [
            0, 
            vertical_size[1], 
            0
            ]
        curve_mesh_rotation = [-np.pi / 2, np.pi / 2, 0]
        self.curve_mesh = Ring(vertical_size[0], outer_radius, mounting_radius, np.pi,
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


class Straight_Spout(ConceptTemplate):
    def __init__(self, num_of_sub_spouts, spout_1_radius, spout_1_thinkness, spout_1_length, spout_1_generatrix_offset, spout_1_rotation, spout_2_radius, spout_2_thinkness, spout_2_length, spout_2_generatrix_offset, spout_2_rotation, spout_3_radius, spout_3_thinkness, spout_3_length, spout_3_generatrix_offset, spout_3_rotation, spout_4_radius, spout_4_thinkness, spout_4_length, spout_4_generatrix_offset, spout_4_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        spout_1_rotation = [x / 180 * np.pi for x in spout_1_rotation]
        spout_2_rotation = [x / 180 * np.pi for x in spout_2_rotation]
        spout_3_rotation = [x / 180 * np.pi for x in spout_3_rotation]
        spout_4_rotation = [x / 180 * np.pi for x in spout_4_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.num_of_sub_spouts = num_of_sub_spouts
        self.spout_1_radius = spout_1_radius
        self.spout_1_thinkness = spout_1_thinkness
        self.spout_1_length = spout_1_length
        self.spout_1_generatrix_offset = spout_1_generatrix_offset
        self.spout_1_rotation = spout_1_rotation
        self.spout_2_radius = spout_2_radius
        self.spout_2_thinkness = spout_2_thinkness
        self.spout_2_length = spout_2_length
        self.spout_2_generatrix_offset = spout_2_generatrix_offset
        self.spout_2_rotation = spout_2_rotation
        self.spout_3_radius = spout_3_radius
        self.spout_3_thinkness = spout_3_thinkness
        self.spout_3_length = spout_3_length
        self.spout_3_generatrix_offset = spout_3_generatrix_offset
        self.spout_3_rotation = spout_3_rotation
        self.spout_4_radius = spout_4_radius
        self.spout_4_thinkness = spout_4_thinkness
        self.spout_4_length = spout_4_length
        self.spout_4_generatrix_offset = spout_4_generatrix_offset
        self.spout_4_rotation = spout_4_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        total_delta_y = -(spout_1_length[0] + spout_1_length[1]) / 4 * np.sin(spout_1_rotation[0])
        total_delta_z = -(spout_1_length[0] + spout_1_length[1]) / 4 * np.cos(spout_1_rotation[0])

        for i in range(num_of_sub_spouts[0]):
            total_delta_y += (locals()['spout_%d_length'%(i+1)][0] + locals()['spout_%d_length'%(i+1)][1]) / 4 * np.sin(locals()['spout_%d_rotation'%(i+1)][0])
            total_delta_z += (locals()['spout_%d_length'%(i+1)][0] + locals()['spout_%d_length'%(i+1)][1]) / 4 * np.cos(locals()['spout_%d_rotation'%(i+1)][0])
            mesh_position = [
                0, 
                total_delta_y, 
                total_delta_z
            ]
            total_delta_y += (locals()['spout_%d_length'%(i+1)][0] + locals()['spout_%d_length'%(i+1)][1]) / 4 * np.sin(locals()['spout_%d_rotation'%(i+1)][0])
            total_delta_z += (locals()['spout_%d_length'%(i+1)][0] + locals()['spout_%d_length'%(i+1)][1]) / 4 * np.cos(locals()['spout_%d_rotation'%(i+1)][0])
            mesh_rotation = [
                np.pi / 2 - locals()['spout_%d_rotation'%(i+1)][0], 
                np.pi / 2, 
                0
            ]
            self.mesh = Ring(height = locals()['spout_%d_length'%(i+1)][0], 
                             outer_top_radius = locals()['spout_%d_radius'%(i+1)][0], 
                             inner_top_radius = locals()['spout_%d_radius'%(i+1)][0] - locals()['spout_%d_thinkness'%(i+1)][0] * 2,
                             outer_bottom_radius = locals()['spout_%d_radius'%(i+1)][1], 
                             inner_bottom_radius = locals()['spout_%d_radius'%(i+1)][1] - locals()['spout_%d_thinkness'%(i+1)][0] * 2,
                             back_height = locals()['spout_%d_length'%(i+1)][1], 
                             generatrix_offset = locals()['spout_%d_generatrix_offset'%(i+1)][0], 
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

        self.semantic = 'Spout'


class Curved_Spout(ConceptTemplate):
    def __init__(self, central_radius, exist_angle, torus_radius, thickness, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        exist_angle = [x / 180 * np.pi for x in exist_angle]
        super().__init__(position, rotation)

        # Record Parameters
        self.central_radius = central_radius
        self.exist_angle = exist_angle
        self.torus_radius = torus_radius
        self.thickness = thickness

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        top_mesh_position = [
            0, 
            central_radius[0], 
            0
            ]
        top_mesh_rotation = [0, 0, -np.pi / 2]
        self.top_mesh = Torus(central_radius[0], torus_radius[0], exist_angle[0], torus_radius[1],
                               position = top_mesh_position,
                               rotation = top_mesh_rotation)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0, 
            central_radius[0] * (1 - np.cos(exist_angle[0])) - central_radius[1] * np.cos(exist_angle[0]), 
            central_radius[0] * np.sin(exist_angle[0]) + central_radius[1] * np.sin(exist_angle[0])
            ]
        bottom_mesh_rotation = [
            -exist_angle[0], 
            0, 
            np.pi / 2
        ]
        self.bottom_mesh = Torus(central_radius[1], torus_radius[1], exist_angle[1], torus_radius[2],
                                 position=bottom_mesh_position,
                                 rotation=bottom_mesh_rotation,
                                 rotation_order = "ZXY")
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Spout'