import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation, adjust_position_from_rotation, list_add
from knowledge_utils import *
import trimesh

class Front_Facing_Roller_Body(ConceptTemplate):
    def __init__(self, outer_size, inner_size, inner_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.outer_size = outer_size
        self.inner_size = inner_size
        self.inner_offset = inner_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            0,
            0,
            -inner_size[1] / 2
        ]
        self.bottom_mesh = Cuboid(outer_size[1], outer_size[0], outer_size[2] - inner_size[1], 
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [
            0,
            0,
            (outer_size[2] - inner_size[1]) / 2
        ]
        self.top_mesh = Box_Cylinder_Ring(outer_size[1], outer_size[0], inner_size[1], inner_size[0], 
                                          inner_cylinder_offset = [inner_offset[0], inner_offset[1]],
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

        self.semantic = 'Body'


class Upright_Roller_Body(ConceptTemplate):
    def __init__(self, outer_size, inner_size, inner_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.outer_size = outer_size
        self.inner_size = inner_size
        self.inner_offset = inner_offset

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
        top_mesh_rotation = [-np.pi / 2, 0, 0]
        self.top_mesh = Box_Cylinder_Ring(outer_size[2], outer_size[0], inner_size[1], inner_size[0], 
                                          inner_cylinder_offset = [inner_offset[0], -inner_offset[1]],
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


class Cuboidal_Body(ConceptTemplate):
    def __init__(self, outer_size, thickness, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.outer_size = outer_size
        self.thickness = thickness

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            0,
            -(outer_size[1] - thickness[2]) / 2,
            0
        ]
        self.bottom_mesh = Cuboid(thickness[2], outer_size[0], outer_size[2], 
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [
            0,
            thickness[2] / 2,
            0
        ]
        self.top_mesh = Rectangular_Ring(outer_size[1] - thickness[2], outer_size[0], outer_size[2], 
                                         outer_size[0] - thickness[1] * 2,
                                         outer_size[2] - thickness[0] * 2,
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

        self.semantic = 'Body'


class Roller_Door(ConceptTemplate):
    def __init__(self, circle_size, middle_size, middle_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.circle_size = circle_size
        self.middle_size = middle_size
        self.middle_offset = middle_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        circle_mesh_position = [
            circle_size[0],
            0,
            circle_size[2] / 2
        ]
        circle_mesh_rotation = [np.pi / 2, 0, 0]
        self.circle_mesh = Ring(circle_size[2], circle_size[0], circle_size[1],
                                position = circle_mesh_position,
                                rotation = circle_mesh_rotation) 
        vertices_list.append(self.circle_mesh.vertices)
        faces_list.append(self.circle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.circle_mesh.vertices)

        middle_mesh_position = [
            circle_size[0],
            0,
            circle_size[2] / 2
        ]
        middle_mesh_rotation = [np.pi / 2, 0, 0]
        self.middle_mesh = Ring(middle_size[1], circle_size[1], middle_size[0],
                                inner_offset = middle_offset,
                                position = middle_mesh_position,
                                rotation = middle_mesh_rotation) 
        vertices_list.append(self.middle_mesh.vertices)
        faces_list.append(self.middle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.middle_mesh.vertices)

        center_mesh_position = [
            circle_size[0],
            middle_offset[0],
            circle_size[2] / 2 + middle_offset[1]
        ]
        center_mesh_rotation = [np.pi / 2, 0, 0]
        self.center_mesh = Cylinder(middle_size[1], middle_size[0],
                                    position = center_mesh_position,
                                    rotation = center_mesh_rotation) 
        vertices_list.append(self.center_mesh.vertices)
        faces_list.append(self.center_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.center_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Door'


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

        self.semantic = 'Door'


class Controller_With_Button(ConceptTemplate):
    def __init__(self, bottom_size, button_1_size, button_1_offset, button_2_size, button_2_offset, button_3_size, button_3_offset, button_4_size, button_4_offset, button_5_size, button_5_offset, button_6_size, button_6_offset, button_7_size, button_7_offset, button_8_size, button_8_offset, num_buttons, position = [0, 0, 0], rotation = [0, 0, 0]):

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