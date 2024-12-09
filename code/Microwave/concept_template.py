import numpy as np
import trimesh
from knowledge_utils import *
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation, adjust_position_from_rotation, list_add

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
        self.vertices = apply_transformation(self.vertices, position, rotation)\
        
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


class Cuboidal_Handle(ConceptTemplate):
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
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="ZXY")

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Trifold_Handle(ConceptTemplate):
    def __init__(self, mounting_size, mounting_seperation, grip_size, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.mounting_size = mounting_size
        self.mounting_seperation = mounting_seperation
        self.grip_size = grip_size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        top_mesh_position = [
            0, 
            mounting_seperation[0] / 2, 
            mounting_size[2] / 2
        ]
        self.top_mesh = Cuboid(mounting_size[1], mounting_size[0], mounting_size[2], 
                               position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0, 
            -mounting_seperation[0] / 2, 
            mounting_size[2] / 2
        ]
        self.bottom_mesh = Cuboid(mounting_size[1], mounting_size[0], mounting_size[2], 
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        vertical_mesh_position = [
            0, 
            0,
            mounting_size[2] + grip_size[2] / 2
        ]
        self.vertical_mesh = Cuboid(grip_size[1], grip_size[0], grip_size[2], 
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


class Trifold_Curve_Handle(ConceptTemplate):
    def __init__(self, mounting_size, mounting_seperation, curve_size, curve_exist_angle, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        curve_exist_angle = [x / 180 * np.pi for x in curve_exist_angle]
        super().__init__(position, rotation)

        # Record Parameters
        self.mounting_size = mounting_size
        self.mounting_seperation = mounting_seperation
        self.curve_size = curve_size
        self.curve_exist_angle = curve_exist_angle

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        top_mesh_position = [
            0, 
            mounting_seperation[0] / 2, 
            mounting_size[2] / 2
        ]
        self.top_mesh = Cuboid(mounting_size[1], mounting_size[0], mounting_size[2], 
                               position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0, 
            -mounting_seperation[0] / 2, 
            mounting_size[2] / 2
        ]
        self.bottom_mesh = Cuboid(mounting_size[1], mounting_size[0], mounting_size[2], 
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        curve_z_offset = mounting_size[2] - np.sqrt(curve_size[1] * curve_size[1] - (mounting_seperation[0] / 2) * (mounting_seperation[0] / 2))
        vertical_mesh_position = [
            0, 
            0,
            curve_z_offset
        ]
        vertical_mesh_rotation = [
            0, 
            -np.pi / 2 + curve_exist_angle[0] / 2,
            np.pi / 2
        ]
        self.vertical_mesh = Ring(curve_size[2], curve_size[0], curve_size[1], curve_exist_angle[0],
                                  position = vertical_mesh_position,
                                  rotation = vertical_mesh_rotation)
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


class Curve_Handle(ConceptTemplate):
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

        vertical_mesh_position = [
            0, 
            0,
            -curve_size[0] * np.cos(curve_exist_angle[0] / 2)
        ]
        vertical_mesh_rotation = [
            0, 
            -np.pi / 2 + curve_exist_angle[0] / 2,
            np.pi / 2
        ]
        self.vertical_mesh = Ring(curve_size[2], curve_size[0], curve_size[1], curve_exist_angle[0],
                                  position = vertical_mesh_position,
                                  rotation = vertical_mesh_rotation)
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


class Controller_With_Button(ConceptTemplate):
    def __init__(self, bottom_size, button_1_size, button_1_offset, button_2_size, button_2_offset, button_3_size, button_3_offset, button_4_size, button_4_offset, num_buttons, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.bottom_size = bottom_size
        self.button_1_size = button_1_size
        self.button_1_offset = button_1_offset
        self.button_2_size = button_2_size
        self.button_2_offset = button_2_offset
        self.button_3_size = button_3_size
        self.button_3_offset = button_3_offset
        self.button_4_size = button_4_size
        self.button_4_offset = button_4_offset
        self.num_buttons = num_buttons

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0,
            0,
            bottom_size[2] / 2
        ]
        self.mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                           position = mesh_position)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        for i in range(num_buttons[0]):
            mesh_position = [
                locals()['button_%d_offset'%(i+1)][0], 
                locals()['button_%d_offset'%(i+1)][1], 
                bottom_size[2] + locals()['button_%d_size'%(i+1)][2] / 2
            ]
            self.mesh = Cuboid(locals()['button_%d_size'%(i+1)][1], locals()['button_%d_size'%(i+1)][0], locals()['button_%d_size'%(i+1)][2], 
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


class Cylindrical_Tray(ConceptTemplate):
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
        self.mesh = Cylinder(size[1], size[0],
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

        self.semantic = 'Tray'