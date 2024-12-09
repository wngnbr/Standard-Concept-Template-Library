import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation, adjust_position_from_rotation, list_add
from knowledge_utils import *
import trimesh

class Cuboidal_Shaft(ConceptTemplate):
    def __init__(self, size, layer_offset, shaft_rotation, up_down_relationship, has_central_shaft, central_shaft_size, central_shaft_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        shaft_rotation = [x / 180 * np.pi for x in shaft_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.layer_offset = layer_offset
        self.shaft_rotation = shaft_rotation
        self.up_down_relationship = up_down_relationship
        self.has_central_shaft = has_central_shaft
        self.central_shaft_size = central_shaft_size
        self.central_shaft_offset = central_shaft_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        if up_down_relationship[0] == 0:
            left_mesh_position = [
                0, 
                size[1] / 2,
                layer_offset[0]
            ]
        elif up_down_relationship[0] == 1:
            left_mesh_position = [
                0, 
                -size[1] / 2,
                layer_offset[0]
            ]
        left_mesh_rotation = [0, shaft_rotation[0], 0]
        self.left_mesh = Cuboid(size[1], size[0], size[2],
                                position = left_mesh_position,
                                rotation = left_mesh_rotation)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        if up_down_relationship[0] == 0:
            right_mesh_position = [
                0, 
                -size[1] / 2,
                -layer_offset[0]
            ]
        elif up_down_relationship[0] == 1:
            right_mesh_position = [
                0, 
                size[1] / 2,
                -layer_offset[0]
            ]
        right_mesh_rotation = [0, -shaft_rotation[0], 0]
        self.right_mesh = Cuboid(size[1], size[0], size[2],
                                 position = right_mesh_position,
                                 rotation = right_mesh_rotation)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        if has_central_shaft[0] == 1:
            mesh_position = [
                central_shaft_offset[0],
                central_shaft_offset[1],
                central_shaft_offset[2]
            ]
            self.mesh = Cylinder(central_shaft_size[1], central_shaft_size[0],
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

        self.semantic = 'Shaft'


class Double_Cuboidal_Shaft(ConceptTemplate):
    def __init__(self, size, front_size, front_offset, layer_offset, shaft_rotation, up_down_relationship, has_central_shaft, central_shaft_size, central_shaft_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        shaft_rotation = [x / 180 * np.pi for x in shaft_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.front_size = front_size
        self.front_offset = front_offset
        self.layer_offset = layer_offset
        self.shaft_rotation = shaft_rotation
        self.up_down_relationship = up_down_relationship
        self.has_central_shaft = has_central_shaft
        self.central_shaft_size = central_shaft_size
        self.central_shaft_offset = central_shaft_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # left
        if up_down_relationship[0] == 0:
            left_mesh_position = [
                0, 
                size[1] / 2,
                layer_offset[0]
            ]
        elif up_down_relationship[0] == 1:
            left_mesh_position = [
                0, 
                -size[1] / 2,
                layer_offset[0]
            ]
        left_mesh_rotation = [0, shaft_rotation[0], 0]
        self.left_mesh = Cuboid(size[1], size[0], size[2],
                                position = left_mesh_position,
                                rotation = left_mesh_rotation)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        # left front
        left_front_mesh_position_1 = [
            (size[0] + front_size[0]) / 2, 
            0,
            front_offset[0]
        ]
        left_front_mesh_rotation = [0, shaft_rotation[0], 0]
        left_front_mesh_position_1 = adjust_position_from_rotation(left_front_mesh_position_1, left_front_mesh_rotation)

        if up_down_relationship[0] == 0:
            left_front_mesh_position_2 = [
                0, 
                size[1] / 2,
                layer_offset[0]
            ]
            left_front_mesh_position = list_add(left_front_mesh_position_1, left_front_mesh_position_2)
        elif up_down_relationship[0] == 1:
            left_front_mesh_position_2 = [
                0, 
                -size[1] / 2,
                layer_offset[0]
            ]
            left_front_mesh_position = list_add(left_front_mesh_position_1, left_front_mesh_position_2)

        self.left_front_mesh = Cuboid(size[1], front_size[0], front_size[1],
                                      position = left_front_mesh_position,
                                      rotation = left_front_mesh_rotation)
        vertices_list.append(self.left_front_mesh.vertices)
        faces_list.append(self.left_front_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_front_mesh.vertices)

        # right
        if up_down_relationship[0] == 0:
            right_mesh_position = [
                0, 
                -size[1] / 2,
                -layer_offset[0]
            ]
        elif up_down_relationship[0] == 1:
            right_mesh_position = [
                0, 
                size[1] / 2,
                -layer_offset[0]
            ]
        right_mesh_rotation = [0, -shaft_rotation[0], 0]
        self.right_mesh = Cuboid(size[1], size[0], size[2],
                                 position = right_mesh_position,
                                 rotation = right_mesh_rotation)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        # right front
        right_front_mesh_position_1 = [
            (size[0] + front_size[0]) / 2, 
            0,
            -front_offset[0]
        ]
        right_front_mesh_rotation = [0, -shaft_rotation[0], 0]
        right_front_mesh_position_1 = adjust_position_from_rotation(right_front_mesh_position_1, right_front_mesh_rotation)

        if up_down_relationship[0] == 0:
            right_front_mesh_position_2 = [
                0, 
                -size[1] / 2,
                -layer_offset[0]
            ]
            right_front_mesh_position = list_add(right_front_mesh_position_1, right_front_mesh_position_2)
        elif up_down_relationship[0] == 1:
            right_front_mesh_position_2 = [
                0, 
                size[1] / 2,
                -layer_offset[0]
            ]
            right_front_mesh_position = list_add(right_front_mesh_position_1, right_front_mesh_position_2)

        self.right_front_mesh = Cuboid(size[1], front_size[0], front_size[1],
                                       position = right_front_mesh_position,
                                       rotation = right_front_mesh_rotation)
        vertices_list.append(self.right_front_mesh.vertices)
        faces_list.append(self.right_front_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_front_mesh.vertices)

        if has_central_shaft[0] == 1:
            mesh_position = [
                central_shaft_offset[0],
                central_shaft_offset[1],
                central_shaft_offset[2]
            ]
            self.mesh = Cylinder(central_shaft_size[1], central_shaft_size[0],
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

        self.semantic = 'Shaft'


class Cylindrical_Shaft(ConceptTemplate):
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

        left_mesh_position = [
            0, 
            size[1] / 4,
            0
        ]
        self.left_mesh = Cylinder(size[1] / 2, size[0], 
                                  position = left_mesh_position)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            0, 
            -size[1] / 4,
            0
        ]
        self.right_mesh = Cylinder(size[1] / 2, size[0], 
                                   position = right_mesh_position)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Shaft'


class Cusp_Blade(ConceptTemplate):
    def __init__(self, root_size, root_z_offset, tip_length, tip_z_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.root_size = root_size
        self.root_z_offset = root_z_offset
        self.tip_length = tip_length
        self.tip_z_offset = tip_z_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        root_mesh_position = [
            -root_size[0] / 2, 
            0,
            0
        ]
        root_mesh_rotation = [np.pi / 2, -np.pi / 2, 0]
        self.root_mesh = Cuboid(root_size[0], root_size[3], root_size[1], 
                                root_size[2], root_size[1], 
                                top_offset = [root_z_offset[0], 0],
                                position = root_mesh_position,
                                rotation = root_mesh_rotation)
        vertices_list.append(self.root_mesh.vertices)
        faces_list.append(self.root_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.root_mesh.vertices)

        tip_mesh_position = [
            -root_size[0] - tip_length[0] / 2, 
            0,
            root_z_offset[0]
        ]
        tip_mesh_rotation = [np.pi / 2, -np.pi / 2, 0]
        self.tip_mesh = Cuboid(tip_length[0], 0, root_size[1], 
                                root_size[3], root_size[1], 
                                top_offset = [tip_z_offset[0], 0],
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

        self.semantic = 'Blade'


class Curved_Blade(ConceptTemplate):
    def __init__(self, root_size, root_z_offset, tip_length, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.root_size = root_size
        self.root_z_offset = root_z_offset
        self.tip_length = tip_length

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        root_mesh_position = [
            -root_size[0] / 2, 
            0,
            0
        ]
        root_mesh_rotation = [np.pi / 2, -np.pi / 2, 0]
        self.root_mesh = Cuboid(root_size[0], root_size[3], root_size[1], 
                                root_size[2], root_size[1], 
                                top_offset = [root_z_offset[0], 0],
                                position = root_mesh_position,
                                rotation = root_mesh_rotation)
        vertices_list.append(self.root_mesh.vertices)
        faces_list.append(self.root_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.root_mesh.vertices)

        tip_mesh_position = [
            -root_size[0], 
            0,
            -root_size[3] / 2 + root_z_offset[0]
        ]
        tip_mesh_rotation = [0, 0, np.pi]
        self.tip_mesh = Cylinder(root_size[1], tip_length[0], tip_length[0],
                                  top_radius_z = root_size[3], bottom_radius_z = root_size[3],
                                  is_quarter = True,
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

        self.semantic = 'Blade'


class Ring_Handle(ConceptTemplate):
    def __init__(self, root_size, root_rotation, arm_radius, arm_thickness, arm_offset, arm_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        root_rotation = [x / 180 * np.pi for x in root_rotation]
        arm_rotation = [x / 180 * np.pi for x in arm_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.root_size = root_size
        self.root_rotation = root_rotation
        self.arm_radius = arm_radius
        self.arm_thickness = arm_thickness
        self.arm_offset = arm_offset
        self.arm_rotation = arm_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # root
        root_mesh_position_1 = [
            root_size[0] / 2, 
            0,
            0
        ]
        root_mesh_rotation_1 = [0, -root_rotation[0], 0]
        root_mesh_position = adjust_position_from_rotation(root_mesh_position_1, root_mesh_rotation_1)

        root_mesh_rotation_2 = [np.pi / 2, np.pi / 2, 0]
        root_mesh_rotation = list_add(root_mesh_rotation_1, root_mesh_rotation_2)

        self.root_mesh = Cuboid(root_size[0], root_size[3], root_size[1],
                                root_size[2], root_size[1],
                                position = root_mesh_position,
                                rotation = root_mesh_rotation)
        vertices_list.append(self.root_mesh.vertices)
        faces_list.append(self.root_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.root_mesh.vertices)

        # arm
        arm_mesh_position_1 = [
            arm_radius[0], 
            0,
            0
        ]
        arm_mesh_rotation = [0, -arm_rotation[0], 0]
        arm_mesh_position_1 = adjust_position_from_rotation(arm_mesh_position_1, arm_mesh_rotation)

        arm_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        arm_mesh_position = list_add(arm_mesh_position_1, arm_mesh_position_2)

        self.arm_mesh = Ring(arm_thickness[0], arm_radius[1], arm_radius[3],
                             x_z_ratio = arm_radius[0] / arm_radius[1],
                             inner_x_z_ratio = arm_radius[2] / arm_radius[3],
                             position = arm_mesh_position,
                             rotation = arm_mesh_rotation)
        vertices_list.append(self.arm_mesh.vertices)
        faces_list.append(self.arm_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.arm_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Half_Ring_Handle(ConceptTemplate):
    def __init__(self, root_size, root_rotation, arm_radius, arm_thickness, arm_horizontal_thickness, arm_offset, arm_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        root_rotation = [x / 180 * np.pi for x in root_rotation]
        arm_rotation = [x / 180 * np.pi for x in arm_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.root_size = root_size
        self.root_rotation = root_rotation
        self.arm_radius = arm_radius
        self.arm_thickness = arm_thickness
        self.arm_horizontal_thickness = arm_horizontal_thickness
        self.arm_offset = arm_offset
        self.arm_rotation = arm_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # root
        root_mesh_position_1 = [
            root_size[0] / 2, 
            0,
            0
        ]
        root_mesh_rotation_1 = [0, -root_rotation[0], 0]
        root_mesh_position = adjust_position_from_rotation(root_mesh_position_1, root_mesh_rotation_1)

        root_mesh_rotation_2 = [np.pi / 2, np.pi / 2, 0]
        root_mesh_rotation = list_add(root_mesh_rotation_1, root_mesh_rotation_2)

        self.root_mesh = Cuboid(root_size[0], root_size[3], root_size[1],
                                root_size[2], root_size[1],
                                position = root_mesh_position,
                                rotation = root_mesh_rotation)
        vertices_list.append(self.root_mesh.vertices)
        faces_list.append(self.root_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.root_mesh.vertices)

        # arm ring
        ring_mesh_position_1 = [
            arm_radius[0], 
            0,
            0
        ]
        ring_mesh_rotation = [0, -arm_rotation[0], 0]
        ring_mesh_position_1 = adjust_position_from_rotation(ring_mesh_position_1, ring_mesh_rotation)

        ring_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        ring_mesh_position = list_add(ring_mesh_position_1, ring_mesh_position_2)

        self.ring_mesh = Ring(arm_thickness[0], arm_radius[1], arm_radius[3], np.pi, 
                             x_z_ratio = arm_radius[0] / arm_radius[1],
                             inner_x_z_ratio = arm_radius[2] / arm_radius[3],
                             position = ring_mesh_position,
                             rotation = ring_mesh_rotation)
        vertices_list.append(self.ring_mesh.vertices)
        faces_list.append(self.ring_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.ring_mesh.vertices)

        # arm horizontal
        horizontal_mesh_position_1 = [
            arm_radius[0], 
            0,
            -arm_horizontal_thickness[0] / 2
        ]
        horizontal_mesh_rotation = [0, -arm_rotation[0], 0]
        horizontal_mesh_position_1 = adjust_position_from_rotation(horizontal_mesh_position_1, horizontal_mesh_rotation)

        horizontal_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        horizontal_mesh_position = list_add(horizontal_mesh_position_1, horizontal_mesh_position_2)

        self.horizontal_mesh = Cuboid(arm_thickness[0], arm_radius[0] * 2, arm_horizontal_thickness[0], 
                                      position = horizontal_mesh_position,
                                      rotation = horizontal_mesh_rotation)
        vertices_list.append(self.horizontal_mesh.vertices)
        faces_list.append(self.horizontal_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.horizontal_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Double_Curved_Handle(ConceptTemplate):
    def __init__(self, root_size, root_rotation, arm_radius, arm_thickness, arm_exist_angle, arm_top_bottom_separation, arm_offset, arm_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        root_rotation = [x / 180 * np.pi for x in root_rotation]
        arm_exist_angle = [x / 180 * np.pi for x in arm_exist_angle]
        arm_rotation = [x / 180 * np.pi for x in arm_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.root_size = root_size
        self.root_rotation = root_rotation
        self.arm_radius = arm_radius
        self.arm_thickness = arm_thickness
        self.arm_exist_angle = arm_exist_angle
        self.arm_top_bottom_separation = arm_top_bottom_separation
        self.arm_offset = arm_offset
        self.arm_rotation = arm_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # root
        root_mesh_position_1 = [
            root_size[0] / 2, 
            0,
            0
        ]
        root_mesh_rotation_1 = [0, -root_rotation[0], 0]
        root_mesh_position = adjust_position_from_rotation(root_mesh_position_1, root_mesh_rotation_1)

        root_mesh_rotation_2 = [np.pi / 2, np.pi / 2, 0]
        root_mesh_rotation = list_add(root_mesh_rotation_1, root_mesh_rotation_2)

        self.root_mesh = Cuboid(root_size[0], root_size[3], root_size[1],
                                root_size[2], root_size[1],
                                position = root_mesh_position,
                                rotation = root_mesh_rotation)
        vertices_list.append(self.root_mesh.vertices)
        faces_list.append(self.root_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.root_mesh.vertices)

        # arm bottom
        bottom_mesh_position_1 = [
            arm_radius[1], 
            0,
            0
        ]
        bottom_mesh_rotation = [0, -arm_rotation[0], 0]
        bottom_mesh_position_1 = adjust_position_from_rotation(bottom_mesh_position_1, bottom_mesh_rotation)

        bottom_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        bottom_mesh_position = list_add(bottom_mesh_position_1, bottom_mesh_position_2)
        bottom_mesh_rotation[1] += -np.pi + arm_exist_angle[1] / 2

        self.bottom_mesh = Ring(arm_thickness[0], arm_radius[1], arm_radius[1] - arm_thickness[1], arm_exist_angle[1], 
                                position = bottom_mesh_position,
                                rotation = bottom_mesh_rotation)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        # arm top
        top_mesh_position_1 = [
            -arm_radius[0] * np.cos(arm_exist_angle[0] / 2) + arm_radius[1] * (1 - np.cos(arm_exist_angle[1] / 2)) + arm_top_bottom_separation[0], 
            0,
            0
        ]
        top_mesh_rotation = [0, -arm_rotation[0], 0]
        top_mesh_position_1 = adjust_position_from_rotation(top_mesh_position_1, top_mesh_rotation)

        top_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        top_mesh_position = list_add(top_mesh_position_1, top_mesh_position_2)
        top_mesh_rotation[1] += arm_exist_angle[0] / 2

        self.top_mesh = Ring(arm_thickness[0], arm_radius[0], arm_radius[0] - arm_thickness[1], arm_exist_angle[0], 
                             position = top_mesh_position,
                             rotation = top_mesh_rotation)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        # arm left
        arm_delta_z = arm_radius[0] * np.sin(arm_exist_angle[0] / 2) - arm_radius[1] * np.sin(arm_exist_angle[1] / 2)
        box_length = np.sqrt(arm_delta_z * arm_delta_z + arm_top_bottom_separation[0] * arm_top_bottom_separation[0])
        box_rotation = np.arctan(arm_delta_z / arm_top_bottom_separation[0])

        left_mesh_position_1 = [
            arm_radius[1] * (1 - np.cos(arm_exist_angle[1] / 2)) + arm_top_bottom_separation[0] / 2 + arm_thickness[1] / 2 * np.sin(box_rotation), 
            0,
            -(arm_radius[0] * np.sin(arm_exist_angle[0] / 2) + arm_radius[1] * np.sin(arm_exist_angle[1] / 2)) / 2 + arm_thickness[1] / 2 * np.cos(box_rotation)
        ]
        left_mesh_rotation = [0, -arm_rotation[0], 0]
        left_mesh_position_1 = adjust_position_from_rotation(left_mesh_position_1, left_mesh_rotation)

        left_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        left_mesh_position = list_add(left_mesh_position_1, left_mesh_position_2)
        left_mesh_rotation[1] += box_rotation

        self.left_mesh = Cuboid(arm_thickness[0], box_length, arm_thickness[1], 
                                position = left_mesh_position,
                                rotation = left_mesh_rotation)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        # arm right
        right_mesh_position_1 = [
            arm_radius[1] * (1 - np.cos(arm_exist_angle[1] / 2)) + arm_top_bottom_separation[0] / 2 + arm_thickness[1] / 2 * np.sin(box_rotation), 
            0,
            (arm_radius[0] * np.sin(arm_exist_angle[0] / 2) + arm_radius[1] * np.sin(arm_exist_angle[1] / 2)) / 2 - arm_thickness[1] / 2 * np.cos(box_rotation)
        ]
        right_mesh_rotation = [0, -arm_rotation[0], 0]
        right_mesh_position_1 = adjust_position_from_rotation(right_mesh_position_1, right_mesh_rotation)

        right_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        right_mesh_position = list_add(right_mesh_position_1, right_mesh_position_2)
        right_mesh_rotation[1] += -box_rotation

        self.right_mesh = Cuboid(arm_thickness[0], box_length, arm_thickness[1], 
                                 position = right_mesh_position,
                                 rotation = right_mesh_rotation)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Triple_Curved_Handle(ConceptTemplate):
    def __init__(self, root_size, root_rotation, arm_radius, arm_thickness, arm_exist_angle, arm_top_bottom_separation, arm_offset, arm_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        root_rotation = [x / 180 * np.pi for x in root_rotation]
        arm_exist_angle = [x / 180 * np.pi for x in arm_exist_angle]
        arm_rotation = [x / 180 * np.pi for x in arm_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.root_size = root_size
        self.root_rotation = root_rotation
        self.arm_radius = arm_radius
        self.arm_thickness = arm_thickness
        self.arm_exist_angle = arm_exist_angle
        self.arm_top_bottom_separation = arm_top_bottom_separation
        self.arm_offset = arm_offset
        self.arm_rotation = arm_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # root
        root_mesh_position_1 = [
            root_size[0] / 2, 
            0,
            0
        ]
        root_mesh_rotation_1 = [0, -root_rotation[0], 0]
        root_mesh_position = adjust_position_from_rotation(root_mesh_position_1, root_mesh_rotation_1)

        root_mesh_rotation_2 = [np.pi / 2, np.pi / 2, 0]
        root_mesh_rotation = list_add(root_mesh_rotation_1, root_mesh_rotation_2)

        self.root_mesh = Cuboid(root_size[0], root_size[3], root_size[1],
                                root_size[2], root_size[1],
                                position = root_mesh_position,
                                rotation = root_mesh_rotation)
        vertices_list.append(self.root_mesh.vertices)
        faces_list.append(self.root_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.root_mesh.vertices)

        # arm bottom
        bottom_mesh_position_1 = [
            arm_radius[0], 
            0,
            0
        ]
        bottom_mesh_rotation = [0, -arm_rotation[0], 0]
        bottom_mesh_position_1 = adjust_position_from_rotation(bottom_mesh_position_1, bottom_mesh_rotation)

        bottom_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        bottom_mesh_position = list_add(bottom_mesh_position_1, bottom_mesh_position_2)
        bottom_mesh_rotation[1] += -np.pi + arm_exist_angle[0] / 2

        self.bottom_mesh = Ring(arm_thickness[0], arm_radius[0], arm_radius[0] - arm_thickness[1], arm_exist_angle[0], 
                                position = bottom_mesh_position,
                                rotation = bottom_mesh_rotation)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        # arm top
        top_mesh_position_1 = [
            -arm_radius[0] * np.cos(arm_exist_angle[0] / 2) + arm_radius[0] * (1 - np.cos(arm_exist_angle[0] / 2)) + arm_top_bottom_separation[0], 
            0,
            0
        ]
        top_mesh_rotation = [0, -arm_rotation[0], 0]
        top_mesh_position_1 = adjust_position_from_rotation(top_mesh_position_1, top_mesh_rotation)

        top_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        top_mesh_position = list_add(top_mesh_position_1, top_mesh_position_2)
        top_mesh_rotation[1] += arm_exist_angle[0] / 2

        self.top_mesh = Ring(arm_thickness[0], arm_radius[0], arm_radius[0] - arm_thickness[1], arm_exist_angle[0], 
                             position = top_mesh_position,
                             rotation = top_mesh_rotation)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        # arm left
        arm_delta_z = arm_radius[0] * np.sin(arm_exist_angle[0] / 2) - arm_radius[0] * np.sin(arm_exist_angle[0] / 2)
        box_length = np.sqrt(arm_delta_z * arm_delta_z + arm_top_bottom_separation[0] * arm_top_bottom_separation[0])
        box_rotation = np.arctan(arm_delta_z / arm_top_bottom_separation[0])

        left_mesh_position_1 = [
            arm_radius[0] * (1 - np.cos(arm_exist_angle[0] / 2)) + arm_top_bottom_separation[0] / 2 + arm_thickness[1] / 2 * np.sin(box_rotation), 
            0,
            -(arm_radius[0] * np.sin(arm_exist_angle[0] / 2) + arm_radius[0] * np.sin(arm_exist_angle[0] / 2)) / 2 + arm_thickness[1] / 2 * np.cos(box_rotation)
        ]
        left_mesh_rotation = [0, -arm_rotation[0], 0]
        left_mesh_position_1 = adjust_position_from_rotation(left_mesh_position_1, left_mesh_rotation)

        left_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        left_mesh_position = list_add(left_mesh_position_1, left_mesh_position_2)
        left_mesh_rotation[1] += box_rotation

        self.left_mesh = Cuboid(arm_thickness[0], box_length, arm_thickness[1], 
                                position = left_mesh_position,
                                rotation = left_mesh_rotation)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        # arm right
        right_exist_angle = np.pi - arm_exist_angle[0]
        right_radius = arm_top_bottom_separation[0] / 2 / np.sin(right_exist_angle / 2)

        right_mesh_position_1 = [
            arm_radius[0] * (1 - np.cos(arm_exist_angle[0] / 2)) + arm_top_bottom_separation[0] / 2, 
            0,
            -right_radius * np.cos(right_exist_angle/ 2) + arm_radius[0] * np.sin(arm_exist_angle[0] / 2)
        ]
        right_mesh_rotation = [0, -arm_rotation[0], 0]
        right_mesh_position_1 = adjust_position_from_rotation(right_mesh_position_1, right_mesh_rotation)

        right_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        right_mesh_position = list_add(right_mesh_position_1, right_mesh_position_2)
        right_mesh_rotation[1] += -np.pi / 2 + right_exist_angle / 2

        self.right_mesh = Ring(arm_thickness[0], right_radius, right_radius - arm_thickness[1], right_exist_angle, 
                               position = right_mesh_position,
                               rotation = right_mesh_rotation)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Cuboidal_Ring_Handle(ConceptTemplate):
    def __init__(self, root_size, root_rotation, arm_length, arm_thickness, arm_top_bottom_offset, arm_top_bottom_rotation, arm_top_bottom_separation, arm_offset, arm_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        root_rotation = [x / 180 * np.pi for x in root_rotation]
        arm_top_bottom_rotation = [x / 180 * np.pi for x in arm_top_bottom_rotation]
        arm_rotation = [x / 180 * np.pi for x in arm_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.root_size = root_size
        self.root_rotation = root_rotation
        self.arm_length = arm_length
        self.arm_thickness = arm_thickness
        self.arm_top_bottom_offset = arm_top_bottom_offset
        self.arm_top_bottom_rotation = arm_top_bottom_rotation
        self.arm_top_bottom_separation = arm_top_bottom_separation
        self.arm_offset = arm_offset
        self.arm_rotation = arm_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # root
        root_mesh_position_1 = [
            root_size[0] / 2, 
            0,
            0
        ]
        root_mesh_rotation_1 = [0, -root_rotation[0], 0]
        root_mesh_position = adjust_position_from_rotation(root_mesh_position_1, root_mesh_rotation_1)

        root_mesh_rotation_2 = [np.pi / 2, np.pi / 2, 0]
        root_mesh_rotation = list_add(root_mesh_rotation_1, root_mesh_rotation_2)

        self.root_mesh = Cuboid(root_size[0], root_size[3], root_size[1],
                                root_size[2], root_size[1],
                                position = root_mesh_position,
                                rotation = root_mesh_rotation)
        vertices_list.append(self.root_mesh.vertices)
        faces_list.append(self.root_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.root_mesh.vertices)

        # arm bottom
        bottom_mesh_position_1 = [
            0, 
            0,
            0
        ]
        bottom_mesh_rotation = [0, -arm_rotation[0], 0]
        bottom_mesh_position_1 = adjust_position_from_rotation(bottom_mesh_position_1, bottom_mesh_rotation)

        bottom_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        bottom_mesh_position = list_add(bottom_mesh_position_1, bottom_mesh_position_2)
        bottom_mesh_rotation[1] += arm_top_bottom_rotation[1]

        self.bottom_mesh = Cuboid(arm_thickness[0], arm_thickness[1], arm_length[1], 
                                  position = bottom_mesh_position,
                                  rotation = bottom_mesh_rotation)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        # arm top
        top_mesh_position_1 = [
            arm_top_bottom_separation[0], 
            0,
            arm_top_bottom_offset[0]
        ]
        top_mesh_rotation = [0, -arm_rotation[0], 0]
        top_mesh_position_1 = adjust_position_from_rotation(top_mesh_position_1, top_mesh_rotation)

        top_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        top_mesh_position = list_add(top_mesh_position_1, top_mesh_position_2)
        top_mesh_rotation[1] += arm_top_bottom_rotation[0]

        self.top_mesh = Cuboid(arm_thickness[0], arm_thickness[1], arm_length[0], 
                               position = top_mesh_position,
                               rotation = top_mesh_rotation)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        # arm left
        left_arm_length_x = arm_top_bottom_separation[0] - arm_length[0] * np.sin(arm_top_bottom_rotation[0]) / 2 + arm_length[1] * np.sin(arm_top_bottom_rotation[1]) / 2
        left_arm_length_z = arm_length[0] * np.cos(arm_top_bottom_rotation[0]) / 2 - arm_length[1] * np.cos(arm_top_bottom_rotation[1]) / 2 - arm_top_bottom_offset[0]
        left_arm_length_t = np.sqrt(left_arm_length_x * left_arm_length_x + left_arm_length_z * left_arm_length_z)
        left_arm_offset_x = (arm_top_bottom_separation[0] - arm_length[0] * np.sin(arm_top_bottom_rotation[0]) / 2 - arm_length[1] * np.sin(arm_top_bottom_rotation[1]) / 2) / 2
        left_arm_offset_z = (arm_top_bottom_offset[0] - arm_length[0] * np.cos(arm_top_bottom_rotation[0]) / 2 - arm_length[1] * np.cos(arm_top_bottom_rotation[1]) / 2) / 2
        left_arm_rotation = np.arctan(left_arm_length_z / left_arm_length_x)

        left_mesh_position_1 = [
            left_arm_offset_x, 
            0,
            left_arm_offset_z
        ]
        left_mesh_rotation = [0, -arm_rotation[0], 0]
        left_mesh_position_1 = adjust_position_from_rotation(left_mesh_position_1, left_mesh_rotation)

        left_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        left_mesh_position = list_add(left_mesh_position_1, left_mesh_position_2)
        left_mesh_rotation[1] += left_arm_rotation

        self.left_mesh = Cuboid(arm_thickness[0], left_arm_length_t, arm_thickness[1], 
                                position = left_mesh_position,
                                rotation = left_mesh_rotation)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        # arm right
        right_arm_length_x = arm_top_bottom_separation[0] + arm_length[0] * np.sin(arm_top_bottom_rotation[0]) / 2 - arm_length[1] * np.sin(arm_top_bottom_rotation[1]) / 2
        right_arm_length_z = arm_length[0] * np.cos(arm_top_bottom_rotation[0]) / 2 - arm_length[1] * np.cos(arm_top_bottom_rotation[1]) / 2 + arm_top_bottom_offset[0]
        right_arm_length_t = np.sqrt(right_arm_length_x * right_arm_length_x + right_arm_length_z * right_arm_length_z)
        right_arm_offset_x = (arm_top_bottom_separation[0] + arm_length[0] * np.sin(arm_top_bottom_rotation[0]) / 2 + arm_length[1] * np.sin(arm_top_bottom_rotation[1]) / 2) / 2
        right_arm_offset_z = (arm_top_bottom_offset[0] + arm_length[0] * np.cos(arm_top_bottom_rotation[0]) / 2 + arm_length[1] * np.cos(arm_top_bottom_rotation[1]) / 2) / 2
        right_arm_rotation = np.arctan(right_arm_length_z / right_arm_length_x)

        right_mesh_position_1 = [
            right_arm_offset_x, 
            0,
            right_arm_offset_z
        ]
        right_mesh_rotation = [0, -arm_rotation[0], 0]
        right_mesh_position_1 = adjust_position_from_rotation(right_mesh_position_1, right_mesh_rotation)

        right_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        right_mesh_position = list_add(right_mesh_position_1, right_mesh_position_2)
        right_mesh_rotation[1] += -right_arm_rotation

        self.right_mesh = Cuboid(arm_thickness[0], right_arm_length_t, arm_thickness[1], 
                                 position = right_mesh_position,
                                 rotation = right_mesh_rotation)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Cuboidal_Handle(ConceptTemplate):
    def __init__(self, root_size, root_rotation, arm_size, arm_offset, arm_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        root_rotation = [x / 180 * np.pi for x in root_rotation]
        arm_rotation = [x / 180 * np.pi for x in arm_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.root_size = root_size
        self.root_rotation = root_rotation
        self.arm_size = arm_size
        self.arm_offset = arm_offset
        self.arm_rotation = arm_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # root
        root_mesh_position_1 = [
            root_size[0] / 2, 
            0,
            0
        ]
        root_mesh_rotation_1 = [0, -root_rotation[0], 0]
        root_mesh_position = adjust_position_from_rotation(root_mesh_position_1, root_mesh_rotation_1)

        root_mesh_rotation_2 = [np.pi / 2, np.pi / 2, 0]
        root_mesh_rotation = list_add(root_mesh_rotation_1, root_mesh_rotation_2)

        self.root_mesh = Cuboid(root_size[0], root_size[3], root_size[1],
                                root_size[2], root_size[1],
                                position = root_mesh_position,
                                rotation = root_mesh_rotation)
        vertices_list.append(self.root_mesh.vertices)
        faces_list.append(self.root_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.root_mesh.vertices)

        # arm
        arm_mesh_position_1 = [
            arm_size[0] / 2, 
            0,
            0
        ]
        arm_mesh_rotation = [0, -arm_rotation[0], 0]
        arm_mesh_position_1 = adjust_position_from_rotation(arm_mesh_position_1, arm_mesh_rotation)

        arm_mesh_position_2 = [
            root_size[0] * np.cos(root_rotation[0]) - arm_offset[0] * np.tan(root_rotation[0]), 
            0,
            root_size[0] * np.sin(root_rotation[0]) + arm_offset[0]
        ]
        arm_mesh_position = list_add(arm_mesh_position_1, arm_mesh_position_2)

        self.arm_mesh = Cuboid(arm_size[1], arm_size[0], arm_size[2],
                               position = arm_mesh_position,
                               rotation = arm_mesh_rotation)
        vertices_list.append(self.arm_mesh.vertices)
        faces_list.append(self.arm_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.arm_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'