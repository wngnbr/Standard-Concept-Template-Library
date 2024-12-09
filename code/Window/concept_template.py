import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh


class Standard_Windowframe(ConceptTemplate):
    def __init__(self, size, inner_size, inner_outer_offset, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.inner_size = inner_size
        self.inner_outer_offset = inner_outer_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        base_mesh_rotation = [np.pi / 2, 0, 0]
        self.base_mesh = Rectangular_Ring(size[2], size[0], size[1],
                                          inner_size[0], inner_size[1],
                                          [inner_outer_offset[0], -inner_outer_offset[1]],
                                          rotation=base_mesh_rotation)
        vertices_list.append(self.base_mesh.vertices)
        faces_list.append(self.base_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.base_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Frame'


class Symmetrical_Window(ConceptTemplate):
    def __init__(self, outside_frame_inner_size, outside_frame_inner_outer_offset, number_of_window, size_0, glass_size_0, glass_offset_0, size_1, glass_size_1, glass_offset_1, size_2, glass_size_2, glass_offset_2, offset_x, offset_z, symmetryOrNot, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.outside_frame_inner_size = outside_frame_inner_size
        self.outside_frame_inner_outer_offset = outside_frame_inner_outer_offset
        self.number_of_window = number_of_window
        self.size_0 = size_0
        self.glass_size_0 = glass_size_0
        self.glass_offset_0 = glass_offset_0
        self.size_1 = size_1
        self.glass_size_1 = glass_size_1
        self.glass_offset_1 = glass_offset_1
        self.size_2 = size_2
        self.glass_size_2 = glass_size_2
        self.glass_offset_2 = glass_offset_2
        self.offset_x = offset_x
        self.offset_z = offset_z
        self.symmetryOrNot = symmetryOrNot

        vertices_list = []
        faces_list = []
        self.meshes = []
        total_num_vertices = 0

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # set z position of three layers
        layer_z_position = [
            -size_0[1] / 2 - size_0[1] / 2 + offset_z[0],
            offset_z[0],
            size_0[1] / 2 + size_1[1] / 2 + offset_z[0],
        ]

        # set asymmetry flag (-1 when glass need flipped)
        asymmetryFlag = 1
        if symmetryOrNot:
            asymmetryFlag = -1
        window_configurations = []

        # window configuration information setting
        window_size = [
            {
                "frame_size": size_0,
                "glass_size": glass_size_0,
                "glass_offset": glass_offset_0,
            },
            {
                "frame_size": size_1,
                "glass_size": glass_size_1,
                "glass_offset": glass_offset_1,
            },
            {
                "frame_size": size_2,
                "glass_size": glass_size_2,
                "glass_offset": glass_offset_2,
            },
        ]
        for size in window_size:
            size["frame_size"] = [
                size["frame_size"][0],
                outside_frame_inner_size[1],
                size["frame_size"][1],
            ]
        if number_of_window[0] > 0:
            if number_of_window[0] == 1:
                window_configurations = [
                    {
                        "window_size": window_size[0],
                        "position": [offset_x[0], 0, layer_z_position[1]],
                        "asymmetryFlag": 1,
                    },
                ]
            elif number_of_window[0] == 2:
                window_configurations = [
                    {
                        "window_size": window_size[0],
                        "position": [offset_x[0], 0, layer_z_position[1]],
                        "asymmetryFlag": 1,
                    },
                    {
                        "window_size": window_size[0],
                        "position": [offset_x[1], 0, layer_z_position[1]],
                        "asymmetryFlag": asymmetryFlag,
                    },
                ]
            elif number_of_window[0] == 3:
                window_configurations = [
                    {
                        "window_size": window_size[0],
                        "position": [offset_x[0], 0, layer_z_position[1]],
                        "asymmetryFlag": 1,
                    },
                    {
                        "window_size": window_size[1],
                        "position": [offset_x[1], 0, layer_z_position[2]],
                        "asymmetryFlag": 1,
                    },
                    {
                        "window_size": window_size[0],
                        "position": [offset_x[2], 0, layer_z_position[1]],
                        "asymmetryFlag": asymmetryFlag,
                    },
                ]
            elif number_of_window[0] == 4:
                window_configurations = [
                    {
                        "window_size": window_size[0],
                        "position": [offset_x[0], 0, layer_z_position[1]],
                        "asymmetryFlag": 1,
                    },
                    {
                        "window_size": window_size[1],
                        "position": [offset_x[1], 0, layer_z_position[2]],
                        "asymmetryFlag": 1,
                    },
                    {
                        "window_size": window_size[1],
                        "position": [offset_x[2], 0, layer_z_position[2]],
                        "asymmetryFlag": asymmetryFlag,
                    },
                    {
                        "window_size": window_size[0],
                        "position": [offset_x[3], 0, layer_z_position[1]],
                        "asymmetryFlag": asymmetryFlag,
                    },
                ]
            elif number_of_window[0] == 5:
                window_configurations = [
                    {
                        "window_size": window_size[2],
                        "position": [offset_x[0], 0, layer_z_position[0]],
                        "asymmetryFlag": 1,
                    },
                    {
                        "window_size": window_size[0],
                        "position": [offset_x[1], 0, layer_z_position[1]],
                        "asymmetryFlag": 1,
                    },
                    {
                        "window_size": window_size[1],
                        "position": [offset_x[2], 0, layer_z_position[2]],
                        "asymmetryFlag": 1,
                    },
                    {
                        "window_size": window_size[0],
                        "position": [offset_x[3], 0, layer_z_position[1]],
                        "asymmetryFlag": asymmetryFlag,
                    },
                    {
                        "window_size": window_size[2],
                        "position": [offset_x[4], 0, layer_z_position[0]],
                        "asymmetryFlag": asymmetryFlag,
                    },
                ]
            elif number_of_window[0] == 6:
                window_configurations = [
                    {
                        "window_size": window_size[2],
                        "position": [offset_x[0], 0, layer_z_position[0]],
                        "asymmetryFlag": 1,
                    },
                    {
                        "window_size": window_size[0],
                        "position": [offset_x[1], 0, layer_z_position[1]],
                        "asymmetryFlag": 1,
                    },
                    {
                        "window_size": window_size[1],
                        "position": [offset_x[2], 0, layer_z_position[2]],
                        "asymmetryFlag": 1,
                    },
                    {
                        "window_size": window_size[1],
                        "position": [offset_x[3], 0, layer_z_position[2]],
                        "asymmetryFlag": asymmetryFlag,
                    },
                    {
                        "window_size": window_size[0],
                        "position": [offset_x[4], 0, layer_z_position[1]],
                        "asymmetryFlag": asymmetryFlag,
                    },
                    {
                        "window_size": window_size[2],
                        "position": [offset_x[5], 0, layer_z_position[0]],
                        "asymmetryFlag": asymmetryFlag,
                    },
                ]

        # flip the glass offset x when symmetry , add offset from the outside frame to window position
        for configuration in window_configurations:
            configuration["window_size"]["glass_offset"][0] = (configuration["window_size"]["glass_offset"][0] * configuration["asymmetryFlag"])
            configuration["position"][1] = (configuration["position"][1] + outside_frame_inner_outer_offset[1])

        self.window_configurations = window_configurations
        
        # meshes definition
        for configuration in window_configurations:
            self.frame_mesh = Rectangular_Ring(
                configuration["window_size"]["frame_size"][2],
                configuration["window_size"]["frame_size"][0],
                configuration["window_size"]["frame_size"][1],
                configuration["window_size"]["glass_size"][0],
                configuration["window_size"]["glass_size"][1],
                [
                    configuration["window_size"]["glass_offset"][0],
                    -configuration["window_size"]["glass_offset"][1],
                ],
                position=configuration["position"],
                rotation=[np.pi / 2, 0, 0],
            )
            vertices_list.append(self.frame_mesh.vertices)
            faces_list.append(self.frame_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.frame_mesh.vertices)

            glass_mesh_position = []
            for i in range(3):
                glass_mesh_position.append(configuration["position"][i] + configuration["window_size"]["glass_offset"][i])

            self.glass_mesh = Cuboid(
                configuration["window_size"]["glass_size"][1],
                configuration["window_size"]["glass_size"][0],
                configuration["window_size"]["glass_size"][2],
                position=glass_mesh_position,
            )
            vertices_list.append(self.glass_mesh.vertices)
            faces_list.append(self.glass_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.glass_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Window'


class Asymmetrical_Window(ConceptTemplate):
    def __init__(self, outside_frame_inner_size, outside_frame_inner_outer_offset, number_of_window, size_0, glass_size_0, glass_offset_0, size_1, glass_size_1, glass_offset_1, size_2, glass_size_2, glass_offset_2, size_3, glass_size_3, glass_offset_3, offset_x, offset_z, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.outside_frame_inner_size = outside_frame_inner_size
        self.outside_frame_inner_outer_offset = outside_frame_inner_outer_offset
        self.number_of_window = number_of_window
        self.size_0 = size_0
        self.glass_size_0 = glass_size_0
        self.glass_offset_0 = glass_offset_0
        self.size_1 = size_1
        self.glass_size_1 = glass_size_1
        self.glass_offset_1 = glass_offset_1
        self.size_2 = size_2
        self.glass_size_2 = glass_size_2
        self.glass_offset_2 = glass_offset_2
        self.size_3 = size_3
        self.glass_size_3 = glass_size_3
        self.glass_offset_3 = glass_offset_3
        self.offset_x = offset_x
        self.offset_z = offset_z

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # set z position of three layers
        layer_z_position = [
            offset_z[0],
            size_0[1] / 2 + size_1[1] / 2 + offset_z[0],
            size_0[1] / 2 + size_1[1] + size_2[1] / 2 + offset_z[0],
            size_0[1] / 2 + size_1[1] + size_2[1] + size_3[1] / 2 + offset_z[0],
        ]
        # window configuration information setting
        window_size = [
            {
                "frame_size": size_0,
                "glass_size": glass_size_0,
                "glass_offset": glass_offset_0,
            },
            {
                "frame_size": size_1,
                "glass_size": glass_size_1,
                "glass_offset": glass_offset_1,
            },
            {
                "frame_size": size_2,
                "glass_size": glass_size_2,
                "glass_offset": glass_offset_2,
            },
            {
                "frame_size": size_3,
                "glass_size": glass_size_3,
                "glass_offset": glass_offset_3,
            },
        ]
        for size in window_size:
            size["frame_size"] = [
                size["frame_size"][0],
                outside_frame_inner_size[1],
                size["frame_size"][1],
            ]
        # set window configurations
        window_configurations = []
        for i in range(number_of_window[0]):
            window_configurations.append(
                {
                    "window_size": window_size[i],
                    "position": [
                        offset_x[i],
                        outside_frame_inner_outer_offset[1],
                        layer_z_position[i],
                    ],
                }
            )
        self.window_configurations = window_configurations

        # meshes definition
        for configuration in window_configurations:
            self.frame_mesh = Rectangular_Ring(
                configuration["window_size"]["frame_size"][2],
                configuration["window_size"]["frame_size"][0],
                configuration["window_size"]["frame_size"][1],
                configuration["window_size"]["glass_size"][0],
                configuration["window_size"]["glass_size"][1],
                [
                    configuration["window_size"]["glass_offset"][0],
                    -configuration["window_size"]["glass_offset"][1],
                ],
                position=configuration["position"],
                rotation=[np.pi / 2, 0, 0],
            )
            vertices_list.append(self.frame_mesh.vertices)
            faces_list.append(self.frame_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.frame_mesh.vertices)

            glass_mesh_position = []
            for i in range(3):
                glass_mesh_position.append(configuration["position"][i] + configuration["window_size"]["glass_offset"][i])

            self.glass_mesh = Cuboid(
                configuration["window_size"]["glass_size"][1],
                configuration["window_size"]["glass_size"][0],
                configuration["window_size"]["glass_size"][2],
                position=glass_mesh_position,
            )
            vertices_list.append(self.glass_mesh.vertices)
            faces_list.append(self.glass_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.glass_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Window'


class VerticalSlid_Window(ConceptTemplate):
    def __init__(self, outside_frame_inner_size, outside_frame_inner_outer_offset, number_of_window, size_0, glass_size_0, glass_offset_0, size_1, glass_size_1, glass_offset_1, offset_y, offset_z, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.outside_frame_inner_size = outside_frame_inner_size
        self.outside_frame_inner_outer_offset = outside_frame_inner_outer_offset
        self.number_of_window = number_of_window
        self.size_0 = size_0
        self.glass_size_0 = glass_size_0
        self.glass_offset_0 = glass_offset_0
        self.size_1 = size_1
        self.glass_size_1 = glass_size_1
        self.glass_offset_1 = glass_offset_1
        self.offset_y = offset_y
        self.offset_z = offset_z

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # set z position of three layers
        layer_z_position = [
            offset_z[0],
            size_0[1] / 2 + size_1[1] / 2 + offset_z[0],
        ]

        # window configuration information setting
        window_size = [
            {
                "frame_size": size_0,
                "glass_size": glass_size_0,
                "glass_offset": glass_offset_0,
            },
            {
                "frame_size": size_1,
                "glass_size": glass_size_1,
                "glass_offset": glass_offset_1,
            },
        ]
        for size in window_size:
            size["frame_size"] = [
                outside_frame_inner_size[0],
                size["frame_size"][0],
                size["frame_size"][1],
            ]

        # set window configurations
        window_configurations = []
        for i in range(number_of_window[0]):
            window_configurations.append(
                {
                    "window_size": window_size[i],
                    "position": [
                        outside_frame_inner_outer_offset[0],
                        offset_y[i],
                        layer_z_position[i],
                    ],
                }
            )
        self.window_configurations = window_configurations

        # meshes definition
        for configuration in window_configurations:
            self.frame_mesh = Rectangular_Ring(
                configuration["window_size"]["frame_size"][2],
                configuration["window_size"]["frame_size"][0],
                configuration["window_size"]["frame_size"][1],
                configuration["window_size"]["glass_size"][0],
                configuration["window_size"]["glass_size"][1],
                [
                    configuration["window_size"]["glass_offset"][0],
                    -configuration["window_size"]["glass_offset"][1],
                ],
                position=configuration["position"],
                rotation=[np.pi / 2, 0, 0],
            )
            vertices_list.append(self.frame_mesh.vertices)
            faces_list.append(self.frame_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.frame_mesh.vertices)

            glass_mesh_position = []
            for i in range(3):
                glass_mesh_position.append(configuration["position"][i] + configuration["window_size"]["glass_offset"][i])

            self.glass_mesh = Cuboid(
                configuration["window_size"]["glass_size"][1],
                configuration["window_size"]["glass_size"][0],
                configuration["window_size"]["glass_size"][2],
                position=glass_mesh_position,
            )
            vertices_list.append(self.glass_mesh.vertices)
            faces_list.append(self.glass_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.glass_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Window'


class Cuboidal_Handle(ConceptTemplate):
    def __init__(self, handle_z_position, window_type, windows_size, num_of_handle, size, offset_x, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.handle_z_position = handle_z_position
        self.window_type = window_type
        self.windows_size = windows_size
        self.num_of_handle = num_of_handle
        self.size = size
        self.offset_x = offset_x

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(num_of_handle[0]):
            position_z = 0
            z_layer_position = handle_z_position[i]
            if window_type == 0:
                if z_layer_position == -1:
                    position_z -= windows_size["size_0"][1] / 2
                if z_layer_position >= 0:
                    position_z += windows_size["size_0"][1] / 2
                if z_layer_position >= 1:
                    position_z += windows_size["size_1"][1]
            elif window_type == 1:
                if z_layer_position >= 0:
                    position_z += windows_size["size_0"][1] / 2
                if z_layer_position >= 1:
                    position_z += windows_size["size_1"][1]
                if z_layer_position >= 2:
                    position_z += windows_size["size_2"][1]
                if z_layer_position == 3:
                    position_z += windows_size["size_3"][1]

            handle_mesh_position = [offset_x[i], 0, position_z + size[2] / 2]
            self.handle_mesh = Cuboid(size[1], size[0], size[2],
                                 position=handle_mesh_position)
            vertices_list.append(self.handle_mesh.vertices)
            faces_list.append(self.handle_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.handle_mesh.vertices)

        for i in range(num_of_handle[1]):
            position_z = 0
            z_layer_position = handle_z_position[i]
            if window_type == 0:
                if z_layer_position == -1:
                    position_z -= (
                        windows_size["size_0"][1] / 2 + windows_size["size_2"][1]
                    )
                if z_layer_position >= 0:
                    position_z += windows_size["size_0"][1] / 2
                if z_layer_position >= 1:
                    position_z += windows_size["size_0"][1]
            elif window_type == 1:
                if z_layer_position >= 0:
                    position_z -= windows_size["size_0"][1] / 2
                if z_layer_position >= 1:
                    position_z += windows_size["size_0"][1]
                if z_layer_position >= 2:
                    position_z += windows_size["size_1"][1]
                if z_layer_position == 3:
                    position_z += windows_size["size_2"][1]

            handle_mesh_position = [offset_x[i + 2], 0, position_z - size[2] / 2]
            self.handle_mesh = Cuboid(size[1], size[0], size[2],
                                 position=handle_mesh_position)
            vertices_list.append(self.handle_mesh.vertices)
            faces_list.append(self.handle_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.handle_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Arched_Handle(ConceptTemplate):
    def __init__(self, handle_z_position, window_type, windows_size, num_of_handle, outer_size, bottom_size, offset_x, seperation, thinner_handle, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.handle_z_position = handle_z_position
        self.window_type = window_type
        self.windows_size = windows_size
        self.num_of_handle = num_of_handle
        self.outer_size = outer_size
        self.bottom_size = bottom_size
        self.offset_x = offset_x
        self.seperation = seperation
        self.thinner_handle = thinner_handle

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        central_angle = (np.arcsin((seperation[0] / 2 + bottom_size[1]) / outer_size[0]) * 2)
        arch_offset_z = bottom_size[2] - np.cos(central_angle / 2) * outer_size[0]

        for i in range(num_of_handle[0]):
            position_z = 0
            z_layer_position = handle_z_position[i]
            if window_type == 0:
                if z_layer_position == -1:
                    position_z -= windows_size["size_0"][1] / 2
                if z_layer_position >= 0:
                    position_z += windows_size["size_0"][1] / 2
                if z_layer_position >= 1:
                    position_z += windows_size["size_1"][1]
            elif window_type == 1:
                if z_layer_position >= 0:
                    position_z += windows_size["size_0"][1] / 2
                if z_layer_position >= 1:
                    position_z += windows_size["size_1"][1]
                if z_layer_position >= 2:
                    position_z += windows_size["size_2"][1]
                if z_layer_position == 3:
                    position_z += windows_size["size_3"][1]

            top_mesh_position = [
                offset_x[i],
                seperation[0] / 2 + bottom_size[1] / 2,
                position_z + bottom_size[2] / 2,
            ]   
            self.top_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                                   position=top_mesh_position)
            vertices_list.append(self.top_mesh.vertices)
            faces_list.append(self.top_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.top_mesh.vertices)

            bottom_mesh_position = [
                offset_x[i],
                -seperation[0] / 2 - bottom_size[1] / 2,
                position_z + bottom_size[2] / 2,
            ]
            self.bottom_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                                      position=bottom_mesh_position)
            vertices_list.append(self.bottom_mesh.vertices)
            faces_list.append(self.bottom_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.bottom_mesh.vertices)

            main_mesh_position = [
                offset_x[i],
                0,
                position_z + +arch_offset_z,
            ]
            main_mesh_rotation = [0, -np.pi / 2 + central_angle / 2, np.pi / 2]
            self.main_mesh = Ring(bottom_size[0], outer_size[0],
                                  seperation[0] / 2 / np.sin(central_angle / 2) + thinner_handle[0],
                                  exist_angle=central_angle,
                                  position=main_mesh_position,
                                  rotation=main_mesh_rotation)
            vertices_list.append(self.main_mesh.vertices)
            faces_list.append(self.main_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.main_mesh.vertices)

        for i in range(num_of_handle[1]):
            position_z = 0
            z_layer_position = handle_z_position[i]
            if window_type == 0:
                if z_layer_position == -1:
                    position_z -= (
                        windows_size["size_0"][1] / 2 + windows_size["size_2"][1]
                    )
                if z_layer_position >= 0:
                    position_z -= windows_size["size_0"][1] / 2
                if z_layer_position >= 1:
                    position_z += windows_size["size_0"][1]
            elif window_type == 1:
                if z_layer_position >= 0:
                    position_z -= windows_size["size_0"][1] / 2
                if z_layer_position >= 1:
                    position_z += windows_size["size_0"][1]
                if z_layer_position >= 2:
                    position_z += windows_size["size_1"][1]
                if z_layer_position == 3:
                    position_z += windows_size["size_2"][1]

            top_mesh_position = [
                offset_x[i+2],
                seperation[0] / 2 + bottom_size[1] / 2,
                position_z - bottom_size[2] / 2,
            ]   
            self.top_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                              position=top_mesh_position)
            vertices_list.append(self.top_mesh.vertices)
            faces_list.append(self.top_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.top_mesh.vertices)

            bottom_mesh_position = [
                offset_x[i+2],
                -seperation[0] / 2 - bottom_size[1] / 2,
                position_z - bottom_size[2] / 2,
            ]
            self.bottom_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                                 position=bottom_mesh_position)
            vertices_list.append(self.bottom_mesh.vertices)
            faces_list.append(self.bottom_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.bottom_mesh.vertices)

            main_mesh_position = [
                offset_x[i+2],
                0,
                position_z - arch_offset_z,
            ]
            main_mesh_rotation = [0, np.pi / 2 + central_angle / 2, np.pi / 2]
            self.main_mesh = Ring(bottom_size[0], outer_size[0],
                                  seperation[0] / 2 / np.sin(central_angle / 2) + thinner_handle[0],
                                  exist_angle=central_angle,
                                  position=main_mesh_position,
                                  rotation=main_mesh_rotation)
            vertices_list.append(self.main_mesh.vertices)
            faces_list.append(self.main_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.main_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class LShaped_Handle(ConceptTemplate):
    def __init__(self, handle_z_position, window_type, windows_size, num_of_handle, size_bottom, size_middle, size_top, offset_middle_y, offset_top_y, offset_x, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.handle_z_position = handle_z_position
        self.window_type = window_type
        self.windows_size = windows_size
        self.num_of_handle = num_of_handle
        self.size_bottom = size_bottom
        self.size_middle = size_middle
        self.size_top = size_top
        self.offset_middle_y = offset_middle_y
        self.offset_top_y = offset_top_y
        self.offset_x = offset_x

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(num_of_handle[0]):
            position_z = 0
            z_layer_position = handle_z_position[i]
            if window_type == 0:
                if z_layer_position == -1:
                    position_z -= windows_size["size_0"][1] / 2
                if z_layer_position >= 0:
                    position_z += windows_size["size_0"][1] / 2
                if z_layer_position >= 1:
                    position_z += windows_size["size_1"][1]
            elif window_type == 1:
                if z_layer_position >= 0:
                    position_z += windows_size["size_0"][1] / 2
                if z_layer_position >= 1:
                    position_z += windows_size["size_1"][1]
                if z_layer_position >= 2:
                    position_z += windows_size["size_2"][1]
                if z_layer_position == 3:
                    position_z += windows_size["size_3"][1]

            bottom_mesh_position = [
                offset_x[i],
                0,
                position_z + size_bottom[2] / 2,
            ]
            self.bottom_mesh = Cuboid(size_bottom[1], size_bottom[0], size_bottom[2],
                                      position=bottom_mesh_position)
            vertices_list.append(self.bottom_mesh.vertices)
            faces_list.append(self.bottom_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.bottom_mesh.vertices)

            middle_mesh_position = [
                offset_x[i],
                offset_middle_y[0],
                position_z + size_bottom[2] + size_middle[2] / 2,
            ]
            self.middle_mesh = Cuboid(size_middle[1], size_middle[0], size_middle[2],
                                      position=middle_mesh_position)
            vertices_list.append(self.middle_mesh.vertices)
            faces_list.append(self.middle_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.middle_mesh.vertices)

            top_mesh_position = [
                offset_x[i],
                offset_middle_y[0] + offset_top_y[0],
                position_z + size_bottom[2] + size_middle[2] + size_top[2] / 2,
            ]
            self.top_mesh = Cuboid(size_top[1], size_top[0], size_top[2],
                                   position=top_mesh_position)
            vertices_list.append(self.top_mesh.vertices)
            faces_list.append(self.top_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.top_mesh.vertices)

        for i in range(num_of_handle[1]):
            position_z = 0
            z_layer_position = handle_z_position[i]
            if window_type == 0:
                if z_layer_position == -1:
                    position_z -= (
                        windows_size["size_0"][1] / 2 + windows_size["size_2"][1]
                    )
                if z_layer_position >= 0:
                    position_z += windows_size["size_0"][1] / 2
                if z_layer_position >= 1:
                    position_z += windows_size["size_0"][1]
            elif window_type == 1:
                if z_layer_position >= 0:
                    position_z -= windows_size["size_0"][1] / 2
                if z_layer_position >= 1:
                    position_z += windows_size["size_0"][1]
                if z_layer_position >= 2:
                    position_z += windows_size["size_1"][1]
                if z_layer_position == 3:
                    position_z += windows_size["size_2"][1]

            bottom_mesh_position = [
                offset_x[i+2],
                0,
                position_z - size_bottom[2] / 2,
            ]  
            self.bottom_mesh = Cuboid(size_bottom[1], size_bottom[0], size_bottom[2],
                                      position=bottom_mesh_position)
            vertices_list.append(self.bottom_mesh.vertices)
            faces_list.append(self.bottom_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.bottom_mesh.vertices)

            middle_mesh_position = [
                offset_x[i+2],
                offset_middle_y[0],
                position_z - size_bottom[2] - size_middle[2] / 2,
            ]
            self.middle_mesh = Cuboid(size_middle[1], size_middle[0], size_middle[2],
                                      position=middle_mesh_position)
            vertices_list.append(self.middle_mesh.vertices)
            faces_list.append(self.middle_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.middle_mesh.vertices)

            top_mesh_position = [
                offset_x[i+2],
                offset_middle_y[0] + offset_top_y[0],
                position_z - size_bottom[2] - size_middle[2] - size_top[2] / 2,
            ]
            self.top_mesh = Cuboid(size_top[1], size_top[0], size_top[2],
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

        self.semantic = 'Handle'