import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation, adjust_position_from_rotation, list_add
from knowledge_utils import *
import trimesh

class Cuboidal_Body(ConceptTemplate):
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

        bottom_mesh_position = [
            0,
            0,
            -(size[2] - thickness[4]) / 2
        ]
        self.bottom_mesh = Cuboid(size[1], size[0], thickness[4], 
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        inner_offset_y = (thickness[1] - thickness[0]) / 2
        inner_offset_x = (thickness[3] - thickness[2]) / 2
        top_mesh_position = [
            0,
            0,
            thickness[4] / 2
        ]
        top_mesh_rotation = [np.pi / 2, 0, 0]
        self.top_mesh = Rectangular_Ring(size[2] - thickness[4], size[0], size[1], 
                                         size[0] - thickness[2] - thickness[3],
                                         size[1] - thickness[0] - thickness[1],
                                         [inner_offset_x, -inner_offset_y],
                                         position = top_mesh_position,
                                         rotation = top_mesh_rotation)
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


class Mutiple_Layer_Body(ConceptTemplate):
    def __init__(self, size, thickness, main_clapboard_size, main_clapboard_offset, num_of_sub_clapboards, sub_clapboard_1_size, sub_clapboard_1_offset, sub_clapboard_2_size, sub_clapboard_2_offset, sub_clapboard_3_size, sub_clapboard_3_offset, sub_clapboard_4_size, sub_clapboard_4_offset, sub_clapboard_5_size, sub_clapboard_5_offset, sub_clapboard_6_size, sub_clapboard_6_offset, sub_clapboard_7_size, sub_clapboard_7_offset, sub_clapboard_8_size, sub_clapboard_8_offset, sub_clapboard_9_size, sub_clapboard_9_offset, sub_clapboard_10_size, sub_clapboard_10_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.thickness = thickness
        self.main_clapboard_size = main_clapboard_size
        self.main_clapboard_offset = main_clapboard_offset
        self.num_of_sub_clapboards = num_of_sub_clapboards
        self.sub_clapboard_1_size = sub_clapboard_1_size
        self.sub_clapboard_1_offset = sub_clapboard_1_offset
        self.sub_clapboard_2_size = sub_clapboard_2_size
        self.sub_clapboard_2_offset = sub_clapboard_2_offset
        self.sub_clapboard_3_size = sub_clapboard_3_size
        self.sub_clapboard_3_offset = sub_clapboard_3_offset
        self.sub_clapboard_4_size = sub_clapboard_4_size
        self.sub_clapboard_4_offset = sub_clapboard_4_offset
        self.sub_clapboard_5_size = sub_clapboard_5_size
        self.sub_clapboard_5_offset = sub_clapboard_5_offset
        self.sub_clapboard_6_size = sub_clapboard_6_size
        self.sub_clapboard_6_offset = sub_clapboard_6_offset
        self.sub_clapboard_7_size = sub_clapboard_7_size
        self.sub_clapboard_7_offset = sub_clapboard_7_offset
        self.sub_clapboard_8_size = sub_clapboard_8_size
        self.sub_clapboard_8_offset = sub_clapboard_8_offset
        self.sub_clapboard_9_size = sub_clapboard_9_size
        self.sub_clapboard_9_offset = sub_clapboard_9_offset
        self.sub_clapboard_10_size = sub_clapboard_10_size
        self.sub_clapboard_10_offset = sub_clapboard_10_offset


        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            0,
            0,
            -(size[2] - thickness[4]) / 2
        ]
        self.bottom_mesh = Cuboid(size[1], size[0], thickness[4], 
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        inner_offset_y = (thickness[1] - thickness[0]) / 2
        inner_offset_x = (thickness[3] - thickness[2]) / 2
        top_mesh_position = [
            0,
            0,
            thickness[4] / 2
        ]
        top_mesh_rotation = [np.pi / 2, 0, 0]
        self.top_mesh = Rectangular_Ring(size[2] - thickness[4], size[0], size[1], 
                                         size[0] - thickness[2] - thickness[3],
                                         size[1] - thickness[0] - thickness[1],
                                         [inner_offset_x, -inner_offset_y],
                                         position = top_mesh_position,
                                         rotation = top_mesh_rotation)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        board_mesh_position = [
            (thickness[3] - thickness[2]) / 2,
            (thickness[1] - thickness[0]) / 2 + main_clapboard_offset[0],
            thickness[4] / 2 - (size[2] - thickness[4] - main_clapboard_size[1]) / 2
        ]
        self.board_mesh = Cuboid(main_clapboard_size[0], size[0] - thickness[2] - thickness[3], main_clapboard_size[1],
                                 position = board_mesh_position)
        vertices_list.append(self.board_mesh.vertices)
        faces_list.append(self.board_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.board_mesh.vertices)

        for i in range(num_of_sub_clapboards[0]):
            board_mesh_position = [
                (thickness[3] - thickness[2]) / 2 + locals()['sub_clapboard_%d_offset'%(i+1)][0],
                (thickness[1] - thickness[0]) / 2 + locals()['sub_clapboard_%d_offset'%(i+1)][1],
                thickness[4] / 2 - (size[2] - thickness[4] - locals()['sub_clapboard_%d_size'%(i+1)][2]) / 2
            ]
            self.board_mesh = Cuboid(locals()['sub_clapboard_%d_size'%(i+1)][1], locals()['sub_clapboard_%d_size'%(i+1)][0], locals()['sub_clapboard_%d_size'%(i+1)][2],
                                     position = board_mesh_position)
            vertices_list.append(self.board_mesh.vertices)
            faces_list.append(self.board_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.board_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Cuboidal_Door(ConceptTemplate):
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
            size[2] / 2
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

        self.semantic = 'Door'


class Behind_Double_Layer_Door(ConceptTemplate):
    def __init__(self, main_size, behind_size, behind_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.main_size = main_size
        self.behind_size = behind_size
        self.behind_offset = behind_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0,
            0,
            main_size[2] / 2
        ]
        self.mesh = Cuboid(main_size[1], main_size[0], main_size[2],
                           position = mesh_position) 
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        behind_mesh_position = [
            behind_offset[0],
            behind_offset[1],
            -behind_size[2] / 2
        ]
        self.behind_mesh = Cuboid(behind_size[1], behind_size[0], behind_size[2],
                                  position = behind_mesh_position) 
        vertices_list.append(self.behind_mesh.vertices)
        faces_list.append(self.behind_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.behind_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Door'


class Front_Double_Layer_Door(ConceptTemplate):
    def __init__(self, main_size, front_size, front_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.main_size = main_size
        self.front_size = front_size
        self.front_offset = front_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0,
            0,
            main_size[2] / 2
        ]
        self.mesh = Cuboid(main_size[1], main_size[0], main_size[2],
                           position = mesh_position) 
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        front_mesh_position = [
            front_offset[0],
            front_offset[1],
            main_size[2] + front_size[2] / 2
        ]
        self.front_mesh = Cuboid(front_size[1], front_size[0], front_size[2],
                                  position = front_mesh_position) 
        vertices_list.append(self.front_mesh.vertices)
        faces_list.append(self.front_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Door'


class Sunken_Door(ConceptTemplate):
    def __init__(self, size, sunken_size, sunken_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.sunken_size = sunken_size
        self.sunken_offset = sunken_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            0,
            0,
            -sunken_size[2] / 2 + size[2] / 2
        ]
        self.bottom_mesh = Cuboid(size[1], size[0], size[2] - sunken_size[2], 
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [
            0,
            0,
            (size[2] - sunken_size[2]) / 2 + size[2] / 2
        ]
        top_mesh_rotation = [np.pi / 2, 0, 0]
        self.top_mesh = Rectangular_Ring(sunken_size[2], size[0], size[1], 
                                         sunken_size[0], sunken_size[1],
                                         [sunken_offset[0], -sunken_offset[1]],
                                         position = top_mesh_position,
                                         rotation = top_mesh_rotation)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Door'


class Cylindrical_Connecter(ConceptTemplate):
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

        self.semantic = 'Connector'


class T_Shaped_Connecter(ConceptTemplate):
    def __init__(self, cylinder_size, lateral_cuboid_size, lateral_cuboid_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.cylinder_size = cylinder_size
        self.lateral_cuboid_size = lateral_cuboid_size
        self.lateral_cuboid_offset = lateral_cuboid_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.mesh = Cylinder(cylinder_size[1], cylinder_size[0])
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        beside_mesh_position = [
            cylinder_size[0] + lateral_cuboid_size[0] / 2,
            lateral_cuboid_offset[0],
            -cylinder_size[0] + lateral_cuboid_size[2] / 2
        ]
        self.beside_mesh = Cuboid(lateral_cuboid_size[1], lateral_cuboid_size[0], lateral_cuboid_size[2], 
                                  position = beside_mesh_position)
        vertices_list.append(self.beside_mesh.vertices)
        faces_list.append(self.beside_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.beside_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Connector'


class Trifold_Handle(ConceptTemplate):
    def __init__(self, bottom_size, bottom_seperation, top_size, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.bottom_size = bottom_size
        self.bottom_seperation = bottom_seperation
        self.top_size = top_size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        top_mesh_position = [
            0, 
            bottom_seperation[0] / 2, 
            bottom_size[2] / 2
        ]
        self.top_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2], 
                               position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0, 
            -bottom_seperation[0] / 2, 
            bottom_size[2] / 2
        ]
        self.bottom_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2], 
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        vertical_mesh_position = [
            0, 
            0,
            bottom_size[2] + top_size[2] / 2
        ]
        self.vertical_mesh = Cuboid(top_size[1], top_size[0], top_size[2], 
                                    position = vertical_mesh_position)
        vertices_list.append(self.vertical_mesh.vertices)
        faces_list.append(self.vertical_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.vertical_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="ZXY")

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Claw_Handle(ConceptTemplate):
    def __init__(self, bottom_size, fork_size, fork_offset, fork_tilt_rotation, num_forks, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        fork_tilt_rotation = [x / 180 * np.pi for x in fork_tilt_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.bottom_size = bottom_size
        self.fork_size = fork_size
        self.fork_offset = fork_offset
        self.fork_tilt_rotation = fork_tilt_rotation
        self.num_forks = num_forks

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            0, 
            0, 
            bottom_size[1] / 2
        ]
        bottom_mesh_rotation = [np.pi / 2, 0, 0]
        self.bottom_mesh = Cylinder(bottom_size[1], bottom_size[0],
                                    position = bottom_mesh_position,
                                    rotation = bottom_mesh_rotation)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        for i in range(num_forks[0]):
            rotation_tmp = np.pi * 2 / num_forks[0] * i
            rotate_length = bottom_size[0] + fork_size[0] / 2 * np.cos(fork_tilt_rotation[0])
            
            mesh_position = [
                rotate_length * np.cos(rotation_tmp), 
                rotate_length * np.sin(rotation_tmp), 
                -fork_size[2] / 2 + bottom_size[1] + fork_size[0] / 2 * np.sin(fork_tilt_rotation[0]) - fork_offset[0]
            ]
            mesh_rotation = [
                0, 
                -fork_tilt_rotation[0], 
                rotation_tmp
            ]
            self.mesh = Cuboid(fork_size[1], fork_size[0], fork_size[2],
                               position = mesh_position,
                               rotation = mesh_rotation)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="ZXY")

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Round_Handle(ConceptTemplate):
    def __init__(self, bottom_size, fork_size, fork_offset, fork_tilt_rotation, circle_size, num_forks, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        fork_tilt_rotation = [x / 180 * np.pi for x in fork_tilt_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.bottom_size = bottom_size
        self.fork_size = fork_size
        self.fork_offset = fork_offset
        self.fork_tilt_rotation = fork_tilt_rotation
        self.circle_size = circle_size
        self.num_forks = num_forks

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            0, 
            0, 
            bottom_size[1] / 2
        ]
        bottom_mesh_rotation = [np.pi / 2, 0, 0]
        self.bottom_mesh = Cylinder(bottom_size[1], bottom_size[0],
                                    position = bottom_mesh_position,
                                    rotation = bottom_mesh_rotation)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        for i in range(num_forks[0]):
            rotation_tmp = np.pi * 2 / num_forks[0] * i
            rotate_length = bottom_size[0] + fork_size[0] / 2 * np.cos(fork_tilt_rotation[0])
            
            mesh_position = [
                rotate_length * np.cos(rotation_tmp), 
                rotate_length * np.sin(rotation_tmp), 
                -fork_size[2] / 2 + bottom_size[1] + fork_size[0] / 2 * np.sin(fork_tilt_rotation[0]) - fork_offset[0]
            ]
            mesh_rotation = [
                0, 
                -fork_tilt_rotation[0], 
                rotation_tmp
            ]
            self.mesh = Cuboid(fork_size[1], fork_size[0], fork_size[2],
                               position = mesh_position,
                               rotation = mesh_rotation)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        outer_radius = circle_size[0]
        inner_radius = bottom_size[0] + fork_size[0] * np.cos(fork_tilt_rotation[0])
        circle_mesh_position = [
            0, 
            0, 
            bottom_size[1] - fork_size[2] / 2 + fork_size[0] * np.sin(fork_tilt_rotation[0]) - fork_offset[0]
        ]
        circle_mesh_rotation = [np.pi / 2, 0, 0]
        self.circle_mesh = Ring(circle_size[1], outer_radius, inner_radius, 
                                position = circle_mesh_position,
                                rotation = circle_mesh_rotation)
        vertices_list.append(self.circle_mesh.vertices)
        faces_list.append(self.circle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.circle_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="ZXY")

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Cylindrical_Dial(ConceptTemplate):
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

        bottom_mesh_position = [
            0, 
            0, 
            bottom_size[2] / 2
        ]
        bottom_mesh_rotation = [-np.pi / 2, 0, 0]
        self.bottom_mesh = Cylinder(bottom_size[2], bottom_size[0], bottom_size[1],
                                    position = bottom_mesh_position,
                                    rotation = bottom_mesh_rotation)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [
            0, 
            0, 
            bottom_size[2] + top_size[1] / 2
        ]
        top_mesh_rotation = [-np.pi / 2, 0, 0]
        self.top_mesh = Cylinder(top_size[1], top_size[0],
                                 position = top_mesh_position,
                                 rotation = top_mesh_rotation)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Dial'


class Regular_Controller(ConceptTemplate):
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

        bottom_mesh_position = [
            0, 
            0, 
            size[2] / 2
        ]
        self.bottom_mesh = Cuboid(size[1], size[0], size[2],
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

        self.semantic = 'Controller'


class Cuboidal_Leg(ConceptTemplate):
    def __init__(self, front_legs_size, rear_legs_size, legs_separation, num_legs, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.front_legs_size = front_legs_size
        self.rear_legs_size = rear_legs_size
        self.legs_separation = legs_separation
        self.num_legs = num_legs

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        if (num_legs[0] == 1):
            mesh_position = [
                0,
                -front_legs_size[1] / 2,
                0
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        elif (num_legs[0] == 2):
            mesh_position = [
                legs_separation[0] / 2, 
                -front_legs_size[1] / 2,
                0
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

            mesh_position = [
                -legs_separation[0] / 2,
                -front_legs_size[1] / 2,
                0
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        elif (num_legs[0] == 3):
            mesh_position = [
                legs_separation[0] / 2,
                -front_legs_size[1] / 2,
                legs_separation[2] / 2
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

            mesh_position = [
                -legs_separation[0] / 2,
                -front_legs_size[1] / 2,
                legs_separation[2] / 2
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

            mesh_position = [
                0,
                -rear_legs_size[1] / 2,
                -legs_separation[2] / 2
            ]
            self.mesh = Cuboid(rear_legs_size[1], rear_legs_size[0], rear_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        elif (num_legs[0] == 4):
            mesh_position = [
                legs_separation[0] / 2,
                -front_legs_size[1] / 2,
                legs_separation[2] / 2
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

            mesh_position = [
                -legs_separation[0] / 2,
                -front_legs_size[1] / 2,
                legs_separation[2] / 2
            ]
            self.mesh = Cuboid(front_legs_size[1], front_legs_size[0], front_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

            mesh_position = [
                legs_separation[1] / 2,
                -rear_legs_size[1] / 2,
                -legs_separation[2] / 2
            ]
            self.mesh = Cuboid(rear_legs_size[1], rear_legs_size[0], rear_legs_size[2],
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

            mesh_position = [
                -legs_separation[1] / 2,
                -rear_legs_size[1] / 2,
                -legs_separation[2] / 2
            ]
            self.mesh = Cuboid(rear_legs_size[1], rear_legs_size[0], rear_legs_size[2],
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

        self.semantic = 'Leg'