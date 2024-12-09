import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation, get_rodrigues_matrix
from knowledge_utils import *
import trimesh


class Storagefurniture_body(ConceptTemplate):
    def __init__(self, size, back_size, left_right_inner_size, base_size, has_lid, lid_size, lid_offset,
                 WHOLE_number_of_layer, WHOLE_layer_sizes, WHOLE_layer_offset, WHOLE_interval_between_layers,
                 storagefurniture_layers_params, additional_layers_params, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        additional_layers_params = [x / 180 * np.pi if ((i > 1) and ((i - 1) % 9 in [6, 7, 8])) else x for i, x in enumerate(additional_layers_params)]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.back_size = back_size
        self.left_right_inner_size = left_right_inner_size
        self.base_size = base_size
        self.has_lid = has_lid
        self.lid_size = lid_size
        self.lid_offset = lid_offset
        self.WHOLE_number_of_layer = WHOLE_number_of_layer
        self.WHOLE_layer_sizes = WHOLE_layer_sizes
        self.WHOLE_layer_offset = WHOLE_layer_offset
        self.WHOLE_interval_between_layers = WHOLE_interval_between_layers
        self.EACH_number_of_layer = [storagefurniture_layers_params[i * 5] for i in range(WHOLE_number_of_layer[0] + 1)]
        self.EACH_layer_sizes = [storagefurniture_layers_params[i * 5 + 1: i * 5 + 3] for i in range(WHOLE_number_of_layer[0] + 1)]
        self.EACH_layer_offset = [storagefurniture_layers_params[i * 5 + 3] for i in range(WHOLE_number_of_layer[0] + 1)]
        self.EACH_interval_between_layers = [storagefurniture_layers_params[i * 5 + 4] for i in range(WHOLE_number_of_layer[0] + 1)]
        self.number_of_additional_layers = additional_layers_params[0]
        self.additional_layers_attributes = additional_layers_params[1:]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        body_mesh_number = 4
        lid_mesh_number = 1 if has_lid[0] == 1 else 0
        WHOLE_layer_mesh_number = WHOLE_number_of_layer[0]
        TOTAL_EACH_layer_mesh_number = sum(self.EACH_number_of_layer)
        TOTAL_mesh_list = [body_mesh_number, lid_mesh_number, WHOLE_layer_mesh_number, TOTAL_EACH_layer_mesh_number]
        record_layer_position = []

        for mesh_idx in range(sum(TOTAL_mesh_list[:3])):
            if mesh_idx < TOTAL_mesh_list[0]:
                if mesh_idx < 2:
                    position_sign = -1 if mesh_idx == 0 else 1
                    mesh_position = [position_sign * (size[0] - left_right_inner_size[0]) / 2, 0, 0]
                    self.mesh = Cuboid(size[1], left_right_inner_size[0], size[2], position=mesh_position)
                elif mesh_idx == 3:
                    mesh_position = [0, -(size[1] - base_size[0]) / 2, 0]
                    self.mesh = Cuboid(base_size[0], size[0] - 2 * left_right_inner_size[0], size[2], position=mesh_position)
                else:
                    mesh_position = [0, 0, (back_size[0] - size[2]) / 2]
                    self.mesh = Cuboid(size[1], size[0], back_size[0], position=mesh_position)
            elif mesh_idx < sum(TOTAL_mesh_list[:2]):
                if has_lid[0] == 1:
                    mesh_position = [lid_offset[0], lid_offset[1] + size[1] / 2 + lid_size[1] / 2, lid_offset[2] + back_size[0] / 2]
                    self.mesh = Cuboid(lid_size[1], lid_size[0], lid_size[2], position=mesh_position)
                else:
                    pass
            else:
                record_layer_position.append(WHOLE_layer_offset[0] + (mesh_idx - sum(TOTAL_mesh_list[:2])) * WHOLE_interval_between_layers[mesh_idx - sum(TOTAL_mesh_list[:2]) - 1])
                mesh_position = [0, size[1] / 2 - record_layer_position[mesh_idx - sum(TOTAL_mesh_list[:2])], 0]
                self.mesh = Cuboid(WHOLE_layer_sizes[0], size[0] - 2 * left_right_inner_size[0], size[2], position=mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        for space_idx in range(WHOLE_layer_mesh_number + 1):
            if WHOLE_layer_mesh_number == 0:
                _height = size[1] - base_size[0]
                _pos = base_size[0] / 2
            elif space_idx == 0:
                _height = WHOLE_layer_offset[0] - WHOLE_layer_sizes[0] / 2
                _pos = size[1] / 2 - WHOLE_layer_offset[0] + _height / 2 + WHOLE_layer_sizes[0] / 2
            elif space_idx == WHOLE_layer_mesh_number:
                _height = size[1] - record_layer_position[space_idx - 1] - WHOLE_layer_sizes[0] / 2 - base_size[0]
                _pos = -size[1] / 2 + _height / 2 + base_size[0]
            else:
                _height = WHOLE_interval_between_layers[space_idx - 1]
                _pos = size[1] / 2 - record_layer_position[space_idx] + _height / 2

            for mesh_idx in range(self.EACH_number_of_layer[space_idx]):
                mesh_position = [-size[0] / 2 + (self.EACH_layer_offset[space_idx] + mesh_idx * self.EACH_interval_between_layers[space_idx]) + left_right_inner_size[0], _pos, 0]
                self.mesh = Cuboid(_height, self.EACH_layer_sizes[space_idx][0], self.EACH_layer_sizes[space_idx][1], position=mesh_position)
                vertices_list.append(self.mesh.vertices)
                faces_list.append(self.mesh.faces + total_num_vertices)
                total_num_vertices += len(self.mesh.vertices)

        for additional_idx in range(self.number_of_additional_layers):
            mesh_position = [self.additional_layers_attributes[9 * additional_idx + 3] - position[0],
                             self.additional_layers_attributes[9 * additional_idx + 4] - position[1],
                             self.additional_layers_attributes[9 * additional_idx + 5] - position[2]]
            mesh_rotation = [self.additional_layers_attributes[9 * additional_idx + 6],
                             self.additional_layers_attributes[9 * additional_idx + 7],
                             self.additional_layers_attributes[9 * additional_idx + 8]]
            self.additional_mesh = Cuboid(self.additional_layers_attributes[9 * additional_idx + 1],
                                          self.additional_layers_attributes[9 * additional_idx],
                                          self.additional_layers_attributes[9 * additional_idx + 2],
                                          position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.additional_mesh.vertices)
            faces_list.append(self.additional_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.additional_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Enclosed_leg(ConceptTemplate):
    def __init__(self, size, inner_sizes, additional_legs_params, position=[0, 0, 0], rotation=[0, 0, 0]):
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        additional_legs_params = [additional_legs_params[0]] + [x / 180 * np.pi if (i - 1) % 9 in [6, 7, 8] else x for i, x in enumerate(additional_legs_params[1:], start=1)]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.inner_sizes = inner_sizes
        self.number_of_additional_legs = additional_legs_params[0]
        self.additional_legs_attributes = additional_legs_params[1:]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for mesh_idx in range(4):
            if size[0] <= 0.01 and size[2] <= 0.01:
                break
            position_sign = -1 if mesh_idx % 2 == 0 else 1
            if mesh_idx < 2:
                mesh_position = [position_sign * (size[0] - inner_sizes[0]) / 2, 0, 0]
                self.mesh = Cuboid(size[1], inner_sizes[0], size[2], position=mesh_position)
            else:
                mesh_position = [0, 0, position_sign * (size[2] - inner_sizes[1]) / 2]
                self.mesh = Cuboid(size[1], size[0], inner_sizes[1], position=mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        for i in range(self.number_of_additional_legs):
            mesh_position = [self.additional_legs_attributes[9 * i + 3],
                             self.additional_legs_attributes[9 * i + 4],
                             self.additional_legs_attributes[9 * i + 5]]
            mesh_rotation = [self.additional_legs_attributes[9 * i + 6],
                             self.additional_legs_attributes[9 * i + 7],
                             self.additional_legs_attributes[9 * i + 8]]
            self.additional_mesh = Cuboid(self.additional_legs_attributes[9 * i + 1],
                                          self.additional_legs_attributes[9 * i],
                                          self.additional_legs_attributes[9 * i + 2],
                                          position=mesh_position, rotation=mesh_rotation)
            vertices_list.append(self.additional_mesh.vertices)
            faces_list.append(self.additional_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.additional_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Leg'


class Regular_door(ConceptTemplate):
    def __init__(self, number_of_door, doors_params, position=[0, 0, 0], rotation=[0, 0, 0]):
        
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        doors_params = [x / 180 * np.pi if i % 12 in [8] else x for i, x in enumerate(doors_params)]
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_door = number_of_door
        self.door_size = [doors_params[i * 12: i * 12 + 3] for i in range(number_of_door[0])]
        self.handle_size = [doors_params[i * 12 + 3: i * 12 + 6] for i in range(number_of_door[0])]
        self.handle_offset = [doors_params[i * 12 + 6: i * 12 + 8] for i in range(number_of_door[0])]
        self.door_rotation = [doors_params[i * 12 + 8] for i in range(number_of_door[0])]
        self.door_offset = [doors_params[i * 12 + 9: i * 12 + 12] for i in range(number_of_door[0])]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for door_idx in range(self.number_of_door[0]):
            for mesh_idx in range(2):
                if mesh_idx == 0:
                    mesh_rotation = [0, self.door_rotation[door_idx], 0]
                    mesh_position = [self.door_offset[door_idx][0],
                                     self.door_offset[door_idx][1],
                                     self.door_offset[door_idx][2]]
                    self.mesh = Cuboid(self.door_size[door_idx][1], self.door_size[door_idx][0], self.door_size[door_idx][2], position=mesh_position, rotation=mesh_rotation)
                else:
                    mesh_rotation = [0, self.door_rotation[door_idx], 0]
                    mesh_position = [
                        self.door_offset[door_idx][0] + self.handle_offset[door_idx][0] * np.cos(self.door_rotation[door_idx]) + self.handle_size[door_idx][2] / 2 * np.sin(
                            self.door_rotation[door_idx]),
                        self.door_offset[door_idx][1] + self.handle_offset[door_idx][1],
                        self.door_offset[door_idx][2] - self.handle_offset[door_idx][0] * np.sin(self.door_rotation[door_idx]) + self.handle_size[door_idx][2] / 2 * np.cos(
                            self.door_rotation[door_idx])]
                    self.mesh = Cuboid(self.handle_size[door_idx][1], self.handle_size[door_idx][0], self.handle_size[door_idx][2], position=mesh_position, rotation=mesh_rotation)

                vertices_list.append(self.mesh.vertices)
                faces_list.append(self.mesh.faces + total_num_vertices)
                total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Door'


class Regular_drawer(ConceptTemplate):
    def __init__(self, number_of_drawer, drawers_params, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_drawer = number_of_drawer
        self.drawer_size = [drawers_params[i * 20: i * 20 + 3] for i in range(number_of_drawer[0])]
        self.bottom_size = [drawers_params[i * 20 + 3] for i in range(number_of_drawer[0])]
        self.front_size = [drawers_params[i * 20 + 4: i * 20 + 7] for i in range(number_of_drawer[0])]
        self.front_offset = [drawers_params[i * 20 + 7] for i in range(number_of_drawer[0])]
        self.left_right_inner_size = [drawers_params[i * 20 + 8] for i in range(number_of_drawer[0])]
        self.rear_front_inner_size = [drawers_params[i * 20 + 9] for i in range(number_of_drawer[0])]
        self.number_of_handle = [drawers_params[i * 20 + 10] for i in range(number_of_drawer[0])]
        self.handle_sizes = [drawers_params[i * 20 + 11: i * 20 + 14] for i in range(number_of_drawer[0])]
        self.handle_offset = [drawers_params[i * 20 + 14: i * 20 + 16] for i in range(number_of_drawer[0])]
        self.handle_separation = [drawers_params[i * 20 + 16] for i in range(number_of_drawer[0])]
        self.drawer_offset = [drawers_params[i * 20 + 17: i * 20 + 20] for i in range(number_of_drawer[0])]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for drawer_idx in range(number_of_drawer[0]):
            for mesh_idx in range(6 + self.number_of_handle[drawer_idx]):
                if mesh_idx < 2:
                    position_sign = -1 if mesh_idx == 0 else 1
                    mesh_position = [position_sign * (
                            self.drawer_size[drawer_idx][0] - self.left_right_inner_size[drawer_idx]) / 2 +
                                     self.drawer_offset[drawer_idx][0],
                                     self.drawer_offset[drawer_idx][1],
                                     self.drawer_offset[drawer_idx][2]]
                    self.mesh = Cuboid(self.drawer_size[drawer_idx][1],
                                       self.left_right_inner_size[drawer_idx],
                                       self.drawer_size[drawer_idx][2],
                                       position=mesh_position)
                elif mesh_idx < 4:
                    position_sign = -1 if mesh_idx == 3 else 1
                    mesh_position = [self.drawer_offset[drawer_idx][0],
                                     self.drawer_offset[drawer_idx][1],
                                     position_sign * (
                                             self.drawer_size[drawer_idx][2] - self.rear_front_inner_size[drawer_idx]) / 2 + self.drawer_offset[drawer_idx][2]]
                    self.mesh = Cuboid(self.drawer_size[drawer_idx][1],
                                       self.drawer_size[drawer_idx][0] - 2 * self.left_right_inner_size[drawer_idx],
                                       self.rear_front_inner_size[drawer_idx],
                                       position=mesh_position)
                elif mesh_idx == 4:
                    mesh_position = [self.drawer_offset[drawer_idx][0],
                                     -self.drawer_size[drawer_idx][1] / 2 + self.drawer_offset[drawer_idx][1] -
                                     self.bottom_size[drawer_idx] / 2,
                                     self.drawer_offset[drawer_idx][2]]
                    self.mesh = Cuboid(self.bottom_size[drawer_idx],
                                       self.drawer_size[drawer_idx][0],
                                       self.drawer_size[drawer_idx][2],
                                       position=mesh_position)
                elif mesh_idx == 5:
                    mesh_position = [self.drawer_offset[drawer_idx][0],
                                     self.drawer_offset[drawer_idx][1] +
                                     self.front_offset[drawer_idx],
                                     self.drawer_offset[drawer_idx][2] + self.drawer_size[drawer_idx][2] / 2 +
                                     self.front_size[drawer_idx][2] / 2]
                    self.mesh = Cuboid(self.front_size[drawer_idx][1],
                                       self.front_size[drawer_idx][0],
                                       self.front_size[drawer_idx][2],
                                       position=mesh_position)
                else:
                    if self.number_of_handle[drawer_idx] == 2:
                        position_sign = 1 if mesh_idx == 6 else -1
                    else:
                        position_sign = 0
                    mesh_position = [self.drawer_offset[drawer_idx][0] + self.handle_offset[drawer_idx][0] +
                                     position_sign * self.handle_separation[drawer_idx] / 2,
                                     self.drawer_offset[drawer_idx][1] + self.handle_offset[drawer_idx][1],
                                     self.drawer_offset[drawer_idx][2] + self.drawer_size[drawer_idx][2] / 2 +
                                     self.front_size[drawer_idx][2] + self.front_size[drawer_idx][2] / 2]
                    self.mesh = Cuboid(self.handle_sizes[drawer_idx][1], self.handle_sizes[drawer_idx][0], self.handle_sizes[drawer_idx][2],
                                       position=mesh_position)
                vertices_list.append(self.mesh.vertices)
                faces_list.append(self.mesh.faces + total_num_vertices)
                total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Drawer'


class Regular_front_panel(ConceptTemplate):
    def __init__(self, number_of_frontPanel, frontPanel_params, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_frontPanel = number_of_frontPanel
        self.frontPanel_size = [frontPanel_params[i * 6: i * 6 + 3] for i in range(number_of_frontPanel[0])]
        self.frontPanel_offset = [frontPanel_params[i * 6 + 3: i * 6 + 6] for i in range(number_of_frontPanel[0])]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for frontPanel_idx in range(number_of_frontPanel[0]):
            mesh_position = [self.frontPanel_offset[frontPanel_idx][0],
                             self.frontPanel_offset[frontPanel_idx][1],
                             self.frontPanel_offset[frontPanel_idx][2]]
            self.mesh = Cuboid(self.frontPanel_size[frontPanel_idx][1],
                               self.frontPanel_size[frontPanel_idx][0],
                               self.frontPanel_size[frontPanel_idx][2],
                               position=mesh_position)

            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Panel'