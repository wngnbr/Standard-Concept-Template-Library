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


class Double_Layer_Body(ConceptTemplate):
    def __init__(self, size, thickness, clapboard_size, clapboard_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.thickness = thickness
        self.clapboard_size = clapboard_size
        self.clapboard_offset = clapboard_offset

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
            (thickness[1] - thickness[0]) / 2 + clapboard_offset[0],
            thickness[4] / 2 - (size[2] - thickness[4] - clapboard_size[1]) / 2
        ]
        self.board_mesh = Cuboid(clapboard_size[0], size[0] - thickness[2] - thickness[3], clapboard_size[1],
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


class Flat_Tray(ConceptTemplate):
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

        self.mesh = Cuboid(size[1], size[0], size[2])
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


class Drawer_Like_Tray(ConceptTemplate):
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

        bottom_mesh_position = [
            0,
            -inner_size[1] / 2,
            0
        ]
        self.bottom_mesh = Cuboid(outer_size[1] - inner_size[1], outer_size[0], outer_size[2], 
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [
            0,
            (outer_size[1] - inner_size[1]) / 2,
            0
        ]
        self.top_mesh = Rectangular_Ring(inner_size[1], outer_size[0], outer_size[2], 
                                         inner_size[0], inner_size[2], 
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

        self.semantic = 'Tray'


class Cuboidal_Baffle(ConceptTemplate):
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

        self.semantic = 'Buffle'


class Controller_With_Button(ConceptTemplate):
    def __init__(self, bottom_size, button_1_size, button_1_offset, button_2_size, button_2_offset, button_3_size, button_3_offset, button_4_size, button_4_offset, button_5_size, button_5_offset, button_6_size, button_6_offset, button_7_size, button_7_offset, button_8_size, button_8_offset, button_9_size, button_9_offset, button_10_size, button_10_offset, num_buttons, position = [0, 0, 0], rotation = [0, 0, 0]):

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
        self.button_5_size = button_5_size
        self.button_5_offset = button_5_offset
        self.button_6_size = button_6_size
        self.button_6_offset = button_6_offset
        self.button_7_size = button_7_size
        self.button_7_offset = button_7_offset
        self.button_8_size = button_8_size
        self.button_8_offset = button_8_offset
        self.button_9_size = button_9_size
        self.button_9_offset = button_9_offset
        self.button_10_size = button_10_size
        self.button_10_offset = button_10_offset
        self.num_buttons = num_buttons

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0,
            0,
            bottom_size[3] / 2
        ]
        self.mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                           bottom_width = bottom_size[3],
                           top_offset = [0, -(bottom_size[3] - bottom_size[2]) / 2],
                           position = mesh_position)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        button_rotation = np.arctan((bottom_size[3] - bottom_size[2]) / bottom_size[1])
        for i in range(num_buttons[0]):
            mesh_position = [
                locals()['button_%d_offset'%(i+1)][0], 
                locals()['button_%d_offset'%(i+1)][1] * np.cos(button_rotation) + locals()['button_%d_size'%(i+1)][2] / 2 * np.sin(button_rotation), 
                (bottom_size[2] + bottom_size[3]) / 2 + locals()['button_%d_size'%(i+1)][2] / 2 * np.cos(button_rotation) - locals()['button_%d_offset'%(i+1)][1] * np.sin(button_rotation)
            ]
            mesh_rotation = [-button_rotation, 0, 0]
            self.mesh = Cuboid(locals()['button_%d_size'%(i+1)][1], locals()['button_%d_size'%(i+1)][0], locals()['button_%d_size'%(i+1)][2], 
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

        self.semantic = 'Button'


class Flat_Top(ConceptTemplate):
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

        self.semantic = 'Body'


class Top_With_Burner(ConceptTemplate):
    def __init__(self, bottom_size, burner_1_size, burner_1_thickness, burner_1_offset, burner_1_central_size, burner_1_central_offset, burner_2_size, burner_2_thickness, burner_2_offset, burner_2_central_size, burner_2_central_offset, burner_3_size, burner_3_thickness, burner_3_offset, burner_3_central_size, burner_3_central_offset, burner_4_size, burner_4_thickness, burner_4_offset, burner_4_central_size, burner_4_central_offset, burner_5_size, burner_5_thickness, burner_5_offset, burner_5_central_size, burner_5_central_offset, burner_6_size, burner_6_thickness, burner_6_offset, burner_6_central_size, burner_6_central_offset, num_burners, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.bottom_size = bottom_size
        self.burner_1_size = burner_1_size
        self.burner_1_thickness = burner_1_thickness
        self.burner_1_offset = burner_1_offset
        self.burner_1_central_size = burner_1_central_size
        self.burner_1_central_offset = burner_1_central_offset
        self.burner_2_size = burner_2_size
        self.burner_2_thickness = burner_2_thickness
        self.burner_2_offset = burner_2_offset
        self.burner_2_central_size = burner_2_central_size
        self.burner_2_central_offset = burner_2_central_offset
        self.burner_3_size = burner_3_size
        self.burner_3_thickness = burner_3_thickness
        self.burner_3_offset = burner_3_offset
        self.burner_3_central_size = burner_3_central_size
        self.burner_3_central_offset = burner_3_central_offset
        self.burner_4_size = burner_4_size
        self.burner_4_thickness = burner_4_thickness
        self.burner_4_offset = burner_4_offset
        self.burner_4_central_size = burner_4_central_size
        self.burner_4_central_offset = burner_4_central_offset
        self.burner_5_size = burner_5_size
        self.burner_5_thickness = burner_5_thickness
        self.burner_5_offset = burner_5_offset
        self.burner_5_central_size = burner_5_central_size
        self.burner_5_central_offset = burner_5_central_offset
        self.burner_6_size = burner_6_size
        self.burner_6_thickness = burner_6_thickness
        self.burner_6_offset = burner_6_offset
        self.burner_6_central_size = burner_6_central_size
        self.burner_6_central_offset = burner_6_central_offset
        self.num_burners = num_burners

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0,
            bottom_size[1] / 2,
            0
        ]
        self.mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                           position = mesh_position)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        for i in range(num_burners[0]):
            mesh_position = [
                locals()['burner_%d_offset'%(i+1)][0], 
                bottom_size[1] + locals()['burner_%d_size'%(i+1)][1] / 2, 
                locals()['burner_%d_offset'%(i+1)][1]
            ]
            self.mesh = Rectangular_Ring(locals()['burner_%d_size'%(i+1)][1], locals()['burner_%d_size'%(i+1)][0], locals()['burner_%d_size'%(i+1)][2], 
                                         locals()['burner_%d_size'%(i+1)][0] - locals()['burner_%d_thickness'%(i+1)][0] * 2,
                                         locals()['burner_%d_size'%(i+1)][2] - locals()['burner_%d_thickness'%(i+1)][0] * 2,
                                         position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

            center_mesh_position = [
                locals()['burner_%d_offset'%(i+1)][0] + locals()['burner_%d_central_offset'%(i+1)][0], 
                bottom_size[1] + locals()['burner_%d_central_size'%(i+1)][1] / 2, 
                locals()['burner_%d_offset'%(i+1)][1] + locals()['burner_%d_central_offset'%(i+1)][1]
            ]
            self.center_mesh = Cylinder(locals()['burner_%d_central_size'%(i+1)][1], locals()['burner_%d_central_size'%(i+1)][0],
                                        position = center_mesh_position)
            vertices_list.append(self.center_mesh.vertices)
            faces_list.append(self.center_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.center_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Burner'


class Multilevel_Leg(ConceptTemplate):
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