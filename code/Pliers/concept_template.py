import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation, adjust_position_from_rotation, list_add
from knowledge_utils import *
import trimesh

class Round_Shaft(ConceptTemplate):
    def __init__(self, size, has_central_shaft, central_shaft_size, central_shaft_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.has_central_shaft = has_central_shaft
        self.central_shaft_size = central_shaft_size
        self.central_shaft_offset = central_shaft_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.mesh = Cylinder(size[1], size[0])
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

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


class Rectangular_Shaft(ConceptTemplate):
    def __init__(self, num_layers, layer_1_size, layer_2_size, layer_2_offset, layer_3_size, layer_3_offset, layer_rotation, has_central_shaft, central_shaft_size, central_shaft_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        layer_rotation = [x / 180 * np.pi for x in layer_rotation]
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.num_layers = num_layers
        self.layer_1_size = layer_1_size
        self.layer_2_size = layer_2_size
        self.layer_2_offset = layer_2_offset
        self.layer_3_size = layer_3_size
        self.layer_3_offset = layer_3_offset
        self.layer_rotation = layer_rotation
        self.has_central_shaft = has_central_shaft
        self.central_shaft_size = central_shaft_size
        self.central_shaft_offset = central_shaft_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        layer_1_mesh_position = [
            0, 
            -layer_1_size[1] / 2,
            0
        ]
        layer_1_mesh_rotation = [0, layer_rotation[0], 0]
        self.layer_1_mesh = Cuboid(layer_1_size[1], layer_1_size[0], layer_1_size[2],
                                   position = layer_1_mesh_position,
                                   rotation = layer_1_mesh_rotation)
        vertices_list.append(self.layer_1_mesh.vertices)
        faces_list.append(self.layer_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.layer_1_mesh.vertices)

        if num_layers[0] >= 2:
            layer_2_mesh_position = [
                layer_2_offset[0], 
                layer_2_size[1] / 2,
                layer_2_offset[1]
            ]
            layer_2_mesh_rotation = [0, layer_rotation[1], 0]
            self.layer_2_mesh = Cuboid(layer_2_size[1], layer_2_size[0], layer_2_size[2],
                                       position = layer_2_mesh_position,
                                       rotation = layer_2_mesh_rotation)
            vertices_list.append(self.layer_2_mesh.vertices)
            faces_list.append(self.layer_2_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.layer_2_mesh.vertices)

        if num_layers[0] >= 3:
            layer_3_mesh_position = [
                layer_3_offset[0], 
                layer_2_size[1] + layer_3_size[1] / 2,
                layer_3_offset[1]
            ]
            layer_3_mesh_rotation = [0, layer_rotation[2], 0]
            self.layer_3_mesh = Cuboid(layer_3_size[1], layer_3_size[0], layer_3_size[2],
                                       position = layer_3_mesh_position,
                                       rotation = layer_3_mesh_rotation)
            vertices_list.append(self.layer_3_mesh.vertices)
            faces_list.append(self.layer_3_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.layer_3_mesh.vertices)

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


class Straight_Handle(ConceptTemplate):
    def __init__(self, front_size, behind_size, handle_separation, handle_rotation, front_behind_offset, left_right_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        handle_rotation = [x / 180 * np.pi for x in handle_rotation]
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.front_size = front_size
        self.behind_size = behind_size
        self.handle_separation = handle_separation
        self.handle_rotation = handle_rotation
        self.front_behind_offset = front_behind_offset
        self.left_right_offset = left_right_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # front_left
        front_1_mesh_position_1 = [
            0, 
            0,
            -front_size[2] / 2
        ]
        front_1_mesh_rotation = [0, -handle_rotation[0], 0]
        front_1_mesh_position_1 = adjust_position_from_rotation(front_1_mesh_position_1, front_1_mesh_rotation)

        front_1_mesh_position_2 = [
            handle_separation[0] / 2, 
            0,
            0
        ]
        front_1_mesh_position = list_add(front_1_mesh_position_1, front_1_mesh_position_2)

        self.front_1_mesh = Cuboid(front_size[1], front_size[0], front_size[2],
                                   position = front_1_mesh_position,
                                   rotation = front_1_mesh_rotation)
        vertices_list.append(self.front_1_mesh.vertices)
        faces_list.append(self.front_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_1_mesh.vertices)

        # behind_left
        behind_1_mesh_position_1 = [
            -behind_size[0] / 2, 
            0,
            -behind_size[2] / 2
        ]
        behind_1_mesh_rotation = [0, -handle_rotation[1], 0]
        behind_1_mesh_position_1 = adjust_position_from_rotation(behind_1_mesh_position_1, behind_1_mesh_rotation)

        behind_1_mesh_position_2 = [
            handle_separation[0] / 2 + front_size[0] / 2 * np.cos(handle_rotation[0]) + front_size[2] * np.sin(handle_rotation[0]) + front_behind_offset[0], 
            front_behind_offset[1],
            -front_size[2] * np.cos(handle_rotation[0]) + front_size[0] / 2 * np.sin(handle_rotation[0])
        ]
        behind_1_mesh_position = list_add(behind_1_mesh_position_1, behind_1_mesh_position_2)

        self.behind_1_mesh = Cuboid(behind_size[1], behind_size[0], behind_size[2],
                                    position = behind_1_mesh_position,
                                    rotation = behind_1_mesh_rotation)
        vertices_list.append(self.behind_1_mesh.vertices)
        faces_list.append(self.behind_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.behind_1_mesh.vertices)

        # front_right
        front_2_mesh_position_1 = [
            0, 
            0,
            -front_size[2] / 2
        ]
        front_2_mesh_rotation = [0, handle_rotation[0], 0]
        front_2_mesh_position_1 = adjust_position_from_rotation(front_2_mesh_position_1, front_2_mesh_rotation)

        front_2_mesh_position_2 = [
            -handle_separation[0] / 2, 
            left_right_offset[0],
            0
        ]
        front_2_mesh_position = list_add(front_2_mesh_position_1, front_2_mesh_position_2)

        self.front_2_mesh = Cuboid(front_size[1], front_size[0], front_size[2],
                                   position = front_2_mesh_position,
                                   rotation = front_2_mesh_rotation)
        vertices_list.append(self.front_2_mesh.vertices)
        faces_list.append(self.front_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_2_mesh.vertices)

        # behind_right
        behind_2_mesh_position_1 = [
            behind_size[0] / 2, 
            0,
            -behind_size[2] / 2
        ]
        behind_2_mesh_rotation = [0, handle_rotation[1], 0]
        behind_2_mesh_position_1 = adjust_position_from_rotation(behind_2_mesh_position_1, behind_2_mesh_rotation)

        behind_2_mesh_position_2 = [
            -handle_separation[0] / 2 - front_size[0] / 2 * np.cos(handle_rotation[0]) - front_size[2] * np.sin(handle_rotation[0]) - front_behind_offset[0], 
            front_behind_offset[1] + left_right_offset[0],
            -front_size[2] * np.cos(handle_rotation[0]) + front_size[0] / 2 * np.sin(handle_rotation[0])
        ]
        behind_2_mesh_position = list_add(behind_2_mesh_position_1, behind_2_mesh_position_2)

        self.behind_2_mesh = Cuboid(behind_size[1], behind_size[0], behind_size[2],
                                    position = behind_2_mesh_position,
                                    rotation = behind_2_mesh_rotation)
        vertices_list.append(self.behind_2_mesh.vertices)
        faces_list.append(self.behind_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.behind_2_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Rear_Curved_Handle(ConceptTemplate):
    def __init__(self, front_size, behind_size, exist_angle, handle_separation, handle_rotation, front_behind_offset, left_right_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        exist_angle = [x / 180 * np.pi for x in exist_angle]
        handle_rotation = [x / 180 * np.pi for x in handle_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.front_size = front_size
        self.behind_size = behind_size
        self.exist_angle = exist_angle
        self.handle_separation = handle_separation
        self.handle_rotation = handle_rotation
        self.front_behind_offset = front_behind_offset
        self.left_right_offset = left_right_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # front_left
        front_1_mesh_position_1 = [
            0, 
            0,
            -front_size[2] / 2
        ]
        front_1_mesh_rotation = [0, -handle_rotation[0], 0]
        front_1_mesh_position_1 = adjust_position_from_rotation(front_1_mesh_position_1, front_1_mesh_rotation)

        front_1_mesh_position_2 = [
            handle_separation[0] / 2, 
            0,
            0
        ]
        front_1_mesh_position = list_add(front_1_mesh_position_1, front_1_mesh_position_2)

        self.front_1_mesh = Cuboid(front_size[1], front_size[0], front_size[2],
                                   position = front_1_mesh_position,
                                   rotation = front_1_mesh_rotation)
        vertices_list.append(self.front_1_mesh.vertices)
        faces_list.append(self.front_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_1_mesh.vertices)

        # behind_left
        behind_1_mesh_rotation_1 = [0, np.pi, np.pi]
        behind_1_mesh_position_1 = [
            -behind_size[1], 
            0,
            0
        ]
        behind_1_mesh_rotation_2 = [0, -handle_rotation[1], 0]
        behind_1_mesh_rotation_2_reverse = [0, handle_rotation[1], 0]
        behind_1_mesh_position_1 = adjust_position_from_rotation(behind_1_mesh_position_1, behind_1_mesh_rotation_2)

        behind_1_mesh_position_2 = [
            handle_separation[0] / 2 - front_size[0] / 2 * np.cos(handle_rotation[0]) + front_size[2] * np.sin(handle_rotation[0]) + front_behind_offset[0], 
            front_behind_offset[1],
            -front_size[2] * np.cos(handle_rotation[0]) + front_size[0] / 2 * np.sin(handle_rotation[0])
        ]
        behind_1_mesh_position = list_add(behind_1_mesh_position_1, behind_1_mesh_position_2)
        behind_1_mesh_rotation = list_add(behind_1_mesh_rotation_1, behind_1_mesh_rotation_2_reverse)

        self.behind_1_mesh = Ring(behind_size[2], behind_size[0], behind_size[1], exist_angle[0],
                                  position = behind_1_mesh_position,
                                  rotation = behind_1_mesh_rotation)
        vertices_list.append(self.behind_1_mesh.vertices)
        faces_list.append(self.behind_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.behind_1_mesh.vertices)

        # front_right
        front_2_mesh_position_1 = [
            0, 
            0,
            -front_size[2] / 2
        ]
        front_2_mesh_rotation = [0, handle_rotation[0], 0]
        front_2_mesh_position_1 = adjust_position_from_rotation(front_2_mesh_position_1, front_2_mesh_rotation)

        front_2_mesh_position_2 = [
            -handle_separation[0] / 2, 
            left_right_offset[0],
            0
        ]
        front_2_mesh_position = list_add(front_2_mesh_position_1, front_2_mesh_position_2)

        self.front_2_mesh = Cuboid(front_size[1], front_size[0], front_size[2],
                                   position = front_2_mesh_position,
                                   rotation = front_2_mesh_rotation)
        vertices_list.append(self.front_2_mesh.vertices)
        faces_list.append(self.front_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_2_mesh.vertices)

        # behind_right
        behind_2_mesh_rotation_1 = [0, np.pi, 0]
        behind_2_mesh_position_1 = [
            behind_size[1], 
            0,
            0
        ]
        behind_2_mesh_rotation_2 = [0, handle_rotation[1], 0]
        behind_2_mesh_position_1 = adjust_position_from_rotation(behind_2_mesh_position_1, behind_2_mesh_rotation_2)

        behind_2_mesh_position_2 = [
            -handle_separation[0] / 2 + front_size[0] / 2 * np.cos(handle_rotation[0]) - front_size[2] * np.sin(handle_rotation[0]) - front_behind_offset[0], 
            front_behind_offset[1],
            -front_size[2] * np.cos(handle_rotation[0]) + front_size[0] / 2 * np.sin(handle_rotation[0])
        ]
        behind_2_mesh_position = list_add(behind_2_mesh_position_1, behind_2_mesh_position_2)
        behind_2_mesh_rotation = list_add(behind_2_mesh_rotation_1, behind_2_mesh_rotation_2)

        self.behind_2_mesh = Ring(behind_size[2], behind_size[0], behind_size[1], exist_angle[0],
                                  position = behind_2_mesh_position,
                                  rotation = behind_2_mesh_rotation)
        vertices_list.append(self.behind_2_mesh.vertices)
        faces_list.append(self.behind_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.behind_2_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Middle_Curved_Handle(ConceptTemplate):
    def __init__(self, front_size, middle_size, exist_angle, behind_size, handle_separation, handle_rotation, front_middle_offset, middle_behind_offset, left_right_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        exist_angle = [x / 180 * np.pi for x in exist_angle]
        handle_rotation = [x / 180 * np.pi for x in handle_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.front_size = front_size
        self.middle_size = middle_size
        self.exist_angle = exist_angle
        self.behind_size = behind_size
        self.handle_separation = handle_separation
        self.handle_rotation = handle_rotation
        self.front_middle_offset = front_middle_offset
        self.middle_behind_offset = middle_behind_offset
        self.left_right_offset = left_right_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # front_left
        front_1_mesh_position_1 = [
            0, 
            0,
            -front_size[2] / 2
        ]
        front_1_mesh_rotation = [0, -handle_rotation[0], 0]
        front_1_mesh_position_1 = adjust_position_from_rotation(front_1_mesh_position_1, front_1_mesh_rotation)

        front_1_mesh_position_2 = [
            handle_separation[0] / 2, 
            0,
            0
        ]
        front_1_mesh_position = list_add(front_1_mesh_position_1, front_1_mesh_position_2)

        self.front_1_mesh = Cuboid(front_size[1], front_size[0], front_size[2],
                                   position = front_1_mesh_position,
                                   rotation = front_1_mesh_rotation)
        vertices_list.append(self.front_1_mesh.vertices)
        faces_list.append(self.front_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_1_mesh.vertices)

        # middle_left
        middle_1_mesh_rotation_1 = [0, np.pi, np.pi]
        middle_1_mesh_position_1 = [
            -middle_size[1], 
            0,
            0
        ]
        middle_1_mesh_rotation_2 = [0, -handle_rotation[1], 0]
        middle_1_mesh_rotation_2_reverse = [0, handle_rotation[1], 0]
        middle_1_mesh_position_1 = adjust_position_from_rotation(middle_1_mesh_position_1, middle_1_mesh_rotation_2)

        middle_1_mesh_position_2 = [
            handle_separation[0] / 2 - front_size[0] / 2 * np.cos(handle_rotation[0]) + front_size[2] * np.sin(handle_rotation[0]) + front_middle_offset[0], 
            front_middle_offset[1],
            -front_size[2] * np.cos(handle_rotation[0]) + front_size[0] / 2 * np.sin(handle_rotation[0])
        ]
        middle_1_mesh_position = list_add(middle_1_mesh_position_1, middle_1_mesh_position_2)
        middle_1_mesh_rotation = list_add(middle_1_mesh_rotation_1, middle_1_mesh_rotation_2_reverse)

        self.middle_1_mesh = Ring(middle_size[2], middle_size[0], middle_size[1], exist_angle[0],
                                  position = middle_1_mesh_position,
                                  rotation = middle_1_mesh_rotation)
        vertices_list.append(self.middle_1_mesh.vertices)
        faces_list.append(self.middle_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.middle_1_mesh.vertices)

        # behind_left
        behind_1_mesh_position_1 = [
            behind_size[0] / 2, 
            0,
            -behind_size[2] / 2
        ]
        behind_1_mesh_rotation = [0, -handle_rotation[2], 0]
        behind_1_mesh_position_1 = adjust_position_from_rotation(behind_1_mesh_position_1, behind_1_mesh_rotation)

        curve_offset_x = middle_size[1] * (1 - np.cos(exist_angle))
        curve_offset_z = middle_size[1] * np.sin(exist_angle)
        middle_offset_x  = curve_offset_x * np.cos(handle_rotation[1]) - curve_offset_z * np.sin(handle_rotation[1])
        middle_offset_z  = curve_offset_x * np.sin(handle_rotation[1]) + curve_offset_z * np.cos(handle_rotation[1])
        behind_1_mesh_position_2 = [
            handle_separation[0] / 2 - front_size[0] / 2 * np.cos(handle_rotation[0]) + front_size[2] * np.sin(handle_rotation[0]) - middle_offset_x + front_middle_offset[0] + middle_behind_offset[0], 
            front_middle_offset[1] + middle_behind_offset[1],
            -front_size[2] * np.cos(handle_rotation[0]) + front_size[0] / 2 * np.sin(handle_rotation[0]) - middle_offset_z
        ]
        behind_1_mesh_position = list_add(behind_1_mesh_position_1, behind_1_mesh_position_2)

        self.behind_1_mesh = Cuboid(behind_size[1], behind_size[0], behind_size[2],
                                    position = behind_1_mesh_position,
                                    rotation = behind_1_mesh_rotation)
        vertices_list.append(self.behind_1_mesh.vertices)
        faces_list.append(self.behind_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.behind_1_mesh.vertices)

        # front_right
        front_2_mesh_position_1 = [
            0, 
            0,
            -front_size[2] / 2
        ]
        front_2_mesh_rotation = [0, handle_rotation[0], 0]
        front_2_mesh_position_1 = adjust_position_from_rotation(front_2_mesh_position_1, front_2_mesh_rotation)

        front_2_mesh_position_2 = [
            -handle_separation[0] / 2, 
            left_right_offset[0],
            0
        ]
        front_2_mesh_position = list_add(front_2_mesh_position_1, front_2_mesh_position_2)

        self.front_2_mesh = Cuboid(front_size[1], front_size[0], front_size[2],
                                   position = front_2_mesh_position,
                                   rotation = front_2_mesh_rotation)
        vertices_list.append(self.front_2_mesh.vertices)
        faces_list.append(self.front_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_2_mesh.vertices)

        # middle_right
        middle_2_mesh_rotation_1 = [0, np.pi, 0]
        middle_2_mesh_position_1 = [
            middle_size[1], 
            0,
            0
        ]
        middle_2_mesh_rotation_2 = [0, handle_rotation[1], 0]
        middle_2_mesh_position_1 = adjust_position_from_rotation(middle_2_mesh_position_1, middle_2_mesh_rotation_2)

        middle_2_mesh_position_2 = [
            -handle_separation[0] / 2 + front_size[0] / 2 * np.cos(handle_rotation[0]) - front_size[2] * np.sin(handle_rotation[0]) - front_middle_offset[0], 
            front_middle_offset[1],
            -front_size[2] * np.cos(handle_rotation[0]) + front_size[0] / 2 * np.sin(handle_rotation[0])
        ]
        middle_2_mesh_position = list_add(middle_2_mesh_position_1, middle_2_mesh_position_2)
        middle_2_mesh_rotation = list_add(middle_2_mesh_rotation_1, middle_2_mesh_rotation_2)

        self.middle_2_mesh = Ring(middle_size[2], middle_size[0], middle_size[1], exist_angle[0],
                                  position = middle_2_mesh_position,
                                  rotation = middle_2_mesh_rotation)
        vertices_list.append(self.middle_2_mesh.vertices)
        faces_list.append(self.middle_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.middle_2_mesh.vertices)

        # behind_right
        behind_2_mesh_position_1 = [
            -behind_size[0] / 2, 
            0,
            -behind_size[2] / 2
        ]
        behind_2_mesh_rotation = [0, handle_rotation[2], 0]
        behind_2_mesh_position_1 = adjust_position_from_rotation(behind_2_mesh_position_1, behind_2_mesh_rotation)

        curve_offset_x = middle_size[1] * (1 - np.cos(exist_angle))
        curve_offset_z = middle_size[1] * np.sin(exist_angle)
        middle_offset_x  = curve_offset_x * np.cos(handle_rotation[1]) - curve_offset_z * np.sin(handle_rotation[1])
        middle_offset_z  = curve_offset_x * np.sin(handle_rotation[1]) + curve_offset_z * np.cos(handle_rotation[1])
        behind_2_mesh_position_2 = [
            -handle_separation[0] / 2 + front_size[0] / 2 * np.cos(handle_rotation[0]) - front_size[2] * np.sin(handle_rotation[0]) + middle_offset_x - front_middle_offset[0] - middle_behind_offset[0], 
            front_middle_offset[1] + middle_behind_offset[1] + left_right_offset[0],
            -front_size[2] * np.cos(handle_rotation[0]) + front_size[0] / 2 * np.sin(handle_rotation[0]) - middle_offset_z
        ]
        behind_2_mesh_position = list_add(behind_2_mesh_position_1, behind_2_mesh_position_2)

        self.behind_2_mesh = Cuboid(behind_size[1], behind_size[0], behind_size[2],
                                    position = behind_2_mesh_position,
                                    rotation = behind_2_mesh_rotation)
        vertices_list.append(self.behind_2_mesh.vertices)
        faces_list.append(self.behind_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.behind_2_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Asymmetric_Straight_Handle(ConceptTemplate):
    def __init__(self, left_front_size, left_behind_size, left_handle_rotation, left_front_behind_offset, right_front_size, right_behind_size, right_handle_rotation, right_front_behind_offset, handle_separation, left_right_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        left_handle_rotation = [x / 180 * np.pi for x in left_handle_rotation]
        right_handle_rotation = [x / 180 * np.pi for x in right_handle_rotation]
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.left_front_size = left_front_size
        self.left_behind_size = left_behind_size
        self.left_handle_rotation = left_handle_rotation
        self.left_front_behind_offset = left_front_behind_offset
        self.right_front_size = right_front_size
        self.right_behind_size = right_behind_size
        self.right_handle_rotation = right_handle_rotation
        self.right_front_behind_offset = right_front_behind_offset
        self.handle_separation = handle_separation
        self.left_right_offset = left_right_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # front_left
        front_1_mesh_position_1 = [
            0, 
            0,
            -left_front_size[2] / 2
        ]
        front_1_mesh_rotation = [0, left_handle_rotation[0], 0]
        front_1_mesh_position_1 = adjust_position_from_rotation(front_1_mesh_position_1, front_1_mesh_rotation)

        front_1_mesh_position_2 = [
            -handle_separation[0] / 2, 
            0,
            0
        ]
        front_1_mesh_position = list_add(front_1_mesh_position_1, front_1_mesh_position_2)

        self.front_1_mesh = Cuboid(left_front_size[1], left_front_size[0], left_front_size[2],
                                   position = front_1_mesh_position,
                                   rotation = front_1_mesh_rotation)
        vertices_list.append(self.front_1_mesh.vertices)
        faces_list.append(self.front_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_1_mesh.vertices)

        # behind_left
        behind_1_mesh_position_1 = [
            left_behind_size[0] / 2, 
            0,
            -left_behind_size[2] / 2
        ]
        behind_1_mesh_rotation = [0, left_handle_rotation[1], 0]
        behind_1_mesh_position_1 = adjust_position_from_rotation(behind_1_mesh_position_1, behind_1_mesh_rotation)

        behind_1_mesh_position_2 = [
            -handle_separation[0] / 2 - left_front_size[0] / 2 * np.cos(left_handle_rotation[0]) - left_front_size[2] * np.sin(left_handle_rotation[0]) - left_front_behind_offset[0], 
            left_front_behind_offset[1],
            -left_front_size[2] * np.cos(left_handle_rotation[0]) + left_front_size[0] / 2 * np.sin(left_handle_rotation[0])
        ]
        behind_1_mesh_position = list_add(behind_1_mesh_position_1, behind_1_mesh_position_2)

        self.behind_1_mesh = Cuboid(left_behind_size[1], left_behind_size[0], left_behind_size[2],
                                    position = behind_1_mesh_position,
                                    rotation = behind_1_mesh_rotation)
        vertices_list.append(self.behind_1_mesh.vertices)
        faces_list.append(self.behind_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.behind_1_mesh.vertices)

        # front_right
        front_2_mesh_position_1 = [
            0, 
            0,
            -right_front_size[2] / 2
        ]
        front_2_mesh_rotation = [0, -right_handle_rotation[0], 0]
        front_2_mesh_position_1 = adjust_position_from_rotation(front_2_mesh_position_1, front_2_mesh_rotation)

        front_2_mesh_position_2 = [
            handle_separation[0] / 2, 
            left_right_offset[0],
            0
        ]
        front_2_mesh_position = list_add(front_2_mesh_position_1, front_2_mesh_position_2)

        self.front_2_mesh = Cuboid(right_front_size[1], right_front_size[0], right_front_size[2],
                                   position = front_2_mesh_position,
                                   rotation = front_2_mesh_rotation)
        vertices_list.append(self.front_2_mesh.vertices)
        faces_list.append(self.front_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_2_mesh.vertices)

        # behind_right
        behind_2_mesh_position_1 = [
            -right_behind_size[0] / 2, 
            0,
            -right_behind_size[2] / 2
        ]
        behind_2_mesh_rotation = [0, -right_handle_rotation[1], 0]
        behind_2_mesh_position_1 = adjust_position_from_rotation(behind_2_mesh_position_1, behind_2_mesh_rotation)

        behind_2_mesh_position_2 = [
            handle_separation[0] / 2 + right_front_size[0] / 2 * np.cos(right_handle_rotation[0]) + right_front_size[2] * np.sin(right_handle_rotation[0]) + right_front_behind_offset[0], 
            right_front_behind_offset[1] + left_right_offset[0],
            -right_front_size[2] * np.cos(right_handle_rotation[0]) + right_front_size[0] / 2 * np.sin(right_handle_rotation[0])
        ]
        behind_2_mesh_position = list_add(behind_2_mesh_position_1, behind_2_mesh_position_2)

        self.behind_2_mesh = Cuboid(right_behind_size[1], right_behind_size[0], right_behind_size[2],
                                    position = behind_2_mesh_position,
                                    rotation = behind_2_mesh_rotation)
        vertices_list.append(self.behind_2_mesh.vertices)
        faces_list.append(self.behind_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.behind_2_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Cusp_Gripper(ConceptTemplate):
    def __init__(self, behind_size, front_size, gripper_separation, gripper_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):
        
        # Process rotation param
        gripper_rotation = [x / 180 * np.pi for x in gripper_rotation]
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.behind_size = behind_size
        self.front_size = front_size
        self.gripper_separation = gripper_separation
        self.gripper_rotation = gripper_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # front_left
        front_1_mesh_rotation_1 = [np.pi / 2, 0, 0]
        front_1_mesh_position_1 = [
            front_size[0] / 2, 
            0,
            behind_size[3] + front_size[3] / 2
        ]
        front_1_mesh_rotation_2 = [0, gripper_rotation[0], 0]
        front_1_mesh_position_1 = adjust_position_from_rotation(front_1_mesh_position_1, front_1_mesh_rotation_2)

        front_1_mesh_position_2 = [
            gripper_separation[0] / 2, 
            0,
            0
        ]
        front_1_mesh_position = list_add(front_1_mesh_position_1, front_1_mesh_position_2)
        front_1_mesh_rotation = list_add(front_1_mesh_rotation_1, front_1_mesh_rotation_2)

        self.front_1_mesh = Cuboid(front_size[3], front_size[1], front_size[2],
                                   front_size[0], front_size[2],
                                   top_offset = [(front_size[1] - front_size[0]) / 2, 0],
                                   position = front_1_mesh_position,
                                   rotation = front_1_mesh_rotation)
        vertices_list.append(self.front_1_mesh.vertices)
        faces_list.append(self.front_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_1_mesh.vertices)

        # behind_left
        behind_1_mesh_rotation_1 = [np.pi / 2, 0, 0]
        behind_1_mesh_position_1 = [
            behind_size[0] / 2, 
            0,
            behind_size[3] / 2
        ]
        behind_1_mesh_rotation_2 = [0, gripper_rotation[0], 0]
        behind_1_mesh_position_1 = adjust_position_from_rotation(behind_1_mesh_position_1, behind_1_mesh_rotation_2)

        behind_1_mesh_position_2 = [
            gripper_separation[0] / 2, 
            0,
            0
        ]
        behind_1_mesh_position = list_add(behind_1_mesh_position_1, behind_1_mesh_position_2)
        behind_1_mesh_rotation = list_add(behind_1_mesh_rotation_1, behind_1_mesh_rotation_2)

        self.behind_1_mesh = Cuboid(behind_size[3], behind_size[1], behind_size[2],
                                    behind_size[0], behind_size[2],
                                    top_offset = [-(behind_size[1] - behind_size[0]) / 2, 0],
                                    position = behind_1_mesh_position,
                                    rotation = behind_1_mesh_rotation)
        vertices_list.append(self.behind_1_mesh.vertices)
        faces_list.append(self.behind_1_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.behind_1_mesh.vertices)

        # front_right
        front_2_mesh_rotation_1 = [np.pi / 2, 0, 0]
        front_2_mesh_position_1 = [
            -front_size[0] / 2, 
            0,
            behind_size[3] + front_size[3] / 2
        ]
        front_2_mesh_rotation_2 = [0, -gripper_rotation[0], 0]
        front_2_mesh_position_1 = adjust_position_from_rotation(front_2_mesh_position_1, front_2_mesh_rotation_2)

        front_2_mesh_position_2 = [
            -gripper_separation[0] / 2, 
            0,
            0
        ]
        front_2_mesh_position = list_add(front_2_mesh_position_1, front_2_mesh_position_2)
        front_2_mesh_rotation = list_add(front_2_mesh_rotation_1, front_2_mesh_rotation_2)

        self.front_2_mesh = Cuboid(front_size[3], front_size[1], front_size[2],
                                   front_size[0], front_size[2],
                                   top_offset = [-(front_size[1] - front_size[0]) / 2, 0],
                                   position = front_2_mesh_position,
                                   rotation = front_2_mesh_rotation)
        vertices_list.append(self.front_2_mesh.vertices)
        faces_list.append(self.front_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_2_mesh.vertices)

        # behind_right
        behind_2_mesh_rotation_1 = [np.pi / 2, 0, 0]
        behind_2_mesh_position_1 = [
            -behind_size[0] / 2, 
            0,
            behind_size[3] / 2
        ]
        behind_2_mesh_rotation_2 = [0, -gripper_rotation[0], 0]
        behind_2_mesh_position_1 = adjust_position_from_rotation(behind_2_mesh_position_1, behind_2_mesh_rotation_2)

        behind_2_mesh_position_2 = [
            -gripper_separation[0] / 2, 
            0,
            0
        ]
        behind_2_mesh_position = list_add(behind_2_mesh_position_1, behind_2_mesh_position_2)
        behind_2_mesh_rotation = list_add(behind_2_mesh_rotation_1, behind_2_mesh_rotation_2)

        self.behind_2_mesh = Cuboid(behind_size[3], behind_size[1], behind_size[2],
                                    behind_size[0], behind_size[2],
                                    top_offset = [(behind_size[1] - behind_size[0]) / 2, 0],
                                    position = behind_2_mesh_position,
                                    rotation = behind_2_mesh_rotation)
        vertices_list.append(self.behind_2_mesh.vertices)
        faces_list.append(self.behind_2_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.behind_2_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Gripper'


class Curved_Gripper(ConceptTemplate):
    def __init__(self, radius, thickness, gripper_separation, gripper_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        gripper_rotation = [x / 180 * np.pi for x in gripper_rotation]
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.radius = radius
        self.thickness = thickness
        self.gripper_separation = gripper_separation
        self.gripper_rotation = gripper_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        left_mesh_position = [
            gripper_separation[0] / 2, 
            0,
            0
        ]
        left_mesh_rotation = [0, gripper_rotation[0], 0]

        self.left_mesh = Cylinder(thickness[0], radius[0], radius[0],
                                  top_radius_z = radius[1], bottom_radius_z = radius[1],
                                  is_quarter = True,
                                  position = left_mesh_position,
                                  rotation = left_mesh_rotation)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            -gripper_separation[0] / 2, 
            0,
            0
        ]
        right_mesh_rotation = [0, -gripper_rotation[0], np.pi]

        self.right_mesh = Cylinder(thickness[0], radius[0], radius[0],
                                   top_radius_z = radius[1], bottom_radius_z = radius[1],
                                   is_quarter = True,
                                   position = right_mesh_position,
                                   rotation = right_mesh_rotation,
                                   rotation_order = 'ZYX')
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Gripper'


class Rectangular_Baffle(ConceptTemplate):
    def __init__(self, size, baffle_separation, baffle_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        baffle_rotation = [x / 180 * np.pi for x in baffle_rotation]
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.baffle_separation = baffle_separation
        self.baffle_rotation = baffle_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        left_mesh_position = [
            baffle_separation[0] / 2, 
            0,
            0
        ]
        left_mesh_rotation = [0, baffle_rotation[0], 0]

        self.left_mesh = Cuboid(size[1], size[0], size[2],
                                position = left_mesh_position,
                                rotation = left_mesh_rotation)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            -baffle_separation[0] / 2, 
            0,
            0
        ]
        right_mesh_rotation = [0, -baffle_rotation[0], 0]

        self.right_mesh = Cuboid(size[1], size[0], size[2],
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

        self.semantic = 'Baffle'


class Curved_Baffle(ConceptTemplate):
    def __init__(self, radius, height, exist_angle, seperation_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        exist_angle = [x / 180 * np.pi for x in exist_angle]
        seperation_rotation = [x / 180 * np.pi for x in seperation_rotation]
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.radius = radius
        self.height = height
        self.exist_angle = exist_angle
        self.seperation_rotation = seperation_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        left_mesh_rotation = [
            0, 
            np.pi / 2 - seperation_rotation[0] / 2, 
            0
        ]

        self.left_mesh = Ring(height[0], radius[0], radius[1], exist_angle[0],
                              rotation = left_mesh_rotation)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_rotation = [
            0, 
            np.pi / 2 - seperation_rotation[0] / 2, 
            np.pi
        ]

        self.right_mesh = Ring(height[0], radius[0], radius[1], exist_angle[0],
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

        self.semantic = 'Baffle'