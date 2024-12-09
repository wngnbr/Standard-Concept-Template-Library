import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh

class Cuboidal_Body(ConceptTemplate):
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

        self.semantic = 'Body'


class Cambered_Body(ConceptTemplate):
    def __init__(self, size, beside_radius_z, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.beside_radius_z = beside_radius_z

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.mesh = Cuboid(size[1], size[0], size[2])
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        left_mesh_position = [
            0, 
            0,
            -size[2] / 2
        ]
        left_mesh_rotation = [0, np.pi, 0]
        self.left_mesh = Cylinder(size[1], size[0] / 2,
                                  top_radius_z = beside_radius_z[0],
                                  bottom_radius_z = beside_radius_z[0],
                                  is_half = True,
                                  position = left_mesh_position,
                                  rotation = left_mesh_rotation)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            0, 
            0,
            size[2] / 2
        ]
        self.right_mesh = Cylinder(size[1], size[0] / 2,
                                   top_radius_z = beside_radius_z[0],
                                   bottom_radius_z = beside_radius_z[0],
                                   is_half = True,
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

        self.semantic = 'Body'


class Double_Layer_Body(ConceptTemplate):
    def __init__(self, main_size, top_size, top_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.main_size = main_size
        self.top_size = top_size
        self.top_offset = top_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.mesh = Cuboid(main_size[1], main_size[0], main_size[2])
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        top_mesh_position = [
            top_offset[0], 
            (main_size[1] + top_size[1]) / 2,
            top_offset[1]
        ]
        self.top_mesh = Cuboid(top_size[1], top_size[0], top_size[2],
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


class Simplified_Wheel(ConceptTemplate):
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

        mesh_rotation = [0, 0, np.pi / 2]
        self.mesh = Cylinder(size[1], size[0], 
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

        self.semantic = 'Wheel'


class Standard_Wheel(ConceptTemplate):
    def __init__(self, middle_size, beside_size, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.middle_size = middle_size
        self.beside_size = beside_size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_rotation = [0, 0, np.pi / 2]
        self.mesh = Cylinder(middle_size[1], middle_size[0], 
                             rotation = mesh_rotation)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        left_mesh_position = [
            (middle_size[1] + beside_size[1]) / 2, 
            0,
            0
        ]
        left_mesh_rotation = [0, 0, np.pi / 2]
        self.left_mesh = Cylinder(beside_size[1], beside_size[0],
                                  position = left_mesh_position,
                                  rotation = left_mesh_rotation)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            -(middle_size[1] + beside_size[1]) / 2, 
            0,
            0
        ]
        right_mesh_rotation = [0, 0, np.pi / 2]
        self.right_mesh = Cylinder(beside_size[1], beside_size[0],
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

        self.semantic = 'Wheel'


class L_Shaped_Button(ConceptTemplate):
    def __init__(self, top_size, bottom_size, top_bottom_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.top_size = top_size
        self.bottom_size = bottom_size
        self.top_bottom_offset = top_bottom_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            0, 
            bottom_size[1] / 2,
            0
        ]
        self.bottom_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [
            0, 
            bottom_size[1] + top_size[1] / 2,
            top_bottom_offset[0]
        ]
        self.top_mesh = Cuboid(top_size[1], top_size[0], top_size[2],
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

        self.semantic = 'Button'


class Double_Cambered_Button(ConceptTemplate):
    def __init__(self, top_size, bottom_size, beside_radius_z, top_bottom_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.top_size = top_size
        self.bottom_size = bottom_size
        self.beside_radius_z = beside_radius_z
        self.top_bottom_offset = top_bottom_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            0, 
            bottom_size[1] / 2,
            0
        ]
        self.bottom_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [
            0, 
            bottom_size[1] + top_size[1] / 2,
            top_bottom_offset[0]
        ]
        self.top_mesh = Cuboid(top_size[1], top_size[0], top_size[2],
                               position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_cambered_mesh_position = [
            0, 
            bottom_size[1] / 2,
            -bottom_size[2] / 2
        ]
        bottom_cambered_mesh_rotation = [0, np.pi, 0]
        self.bottom_cambered_mesh = Cylinder(bottom_size[1], bottom_size[0] / 2,
                                             top_radius_z = beside_radius_z[0],
                                             bottom_radius_z = beside_radius_z[0],
                                             is_half = True,
                                             position = bottom_cambered_mesh_position,
                                             rotation = bottom_cambered_mesh_rotation)
        vertices_list.append(self.bottom_cambered_mesh.vertices)
        faces_list.append(self.bottom_cambered_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_cambered_mesh.vertices)

        top_cambered_mesh_position = [
            0, 
            bottom_size[1] + top_size[1] / 2,
            -top_size[2] / 2 + top_bottom_offset[0]
        ]
        top_cambered_mesh_rotation = [0, np.pi, 0]
        self.top_cambered_mesh = Cylinder(top_size[1], top_size[0] / 2,
                                          top_radius_z = beside_radius_z[1],
                                          bottom_radius_z = beside_radius_z[1],
                                          is_half = True,
                                          position = top_cambered_mesh_position,
                                          rotation = top_cambered_mesh_rotation)
        vertices_list.append(self.top_cambered_mesh.vertices)
        faces_list.append(self.top_cambered_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_cambered_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Button'


class Cuboidal_Nozzle(ConceptTemplate):
    def __init__(self, size, thickness, top_length, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.thickness = thickness
        self.top_length = top_length

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        left_mesh_position = [
            size[0] / 2 - thickness[0] / 2, 
            size[1] / 2,
            0
        ]
        self.left_mesh = Cuboid(size[1], thickness[0], size[2],
                                position = left_mesh_position)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            -size[0] / 2 + thickness[0] / 2, 
            size[1] / 2,
            0
        ]
        self.right_mesh = Cuboid(size[1], thickness[0], size[2],
                                 position = right_mesh_position)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        front_mesh_position = [
            0, 
            size[1] / 2,
            size[2] / 2 - thickness[0] / 2
        ]
        self.front_mesh = Cuboid(size[1], size[0], thickness[0],
                                 position = front_mesh_position)
        vertices_list.append(self.front_mesh.vertices)
        faces_list.append(self.front_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_mesh.vertices)

        top_mesh_position = [
            0, 
            size[1] - thickness[0] / 2,
            (size[2] - top_length[0]) / 2
        ]
        self.top_mesh = Cuboid(thickness[0], size[0], top_length[0],
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

        self.semantic = 'Nozzle'


class Cambered_Nozzle(ConceptTemplate):
    def __init__(self, size, beside_radius_z, thickness, top_length, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.beside_radius_z = beside_radius_z
        self.thickness = thickness
        self.top_length = top_length

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        left_mesh_position = [
            size[0] / 2 - thickness[0] / 2, 
            size[1] / 2,
            0
        ]
        self.left_mesh = Cuboid(size[1], thickness[0], size[2],
                                position = left_mesh_position)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            -size[0] / 2 + thickness[0] / 2, 
            size[1] / 2,
            0
        ]
        self.right_mesh = Cuboid(size[1], thickness[0], size[2],
                                 position = right_mesh_position)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        front_mesh_position = [
            0, 
            size[1] / 2,
            size[2] / 2
        ]
        self.front_mesh = Cylinder(size[1], size[0] / 2,
                                   top_radius_z = beside_radius_z[0],
                                   bottom_radius_z = beside_radius_z[0],
                                   is_half = True,
                                   position = front_mesh_position)
        vertices_list.append(self.front_mesh.vertices)
        faces_list.append(self.front_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_mesh.vertices)

        top_mesh_position = [
            0, 
            size[1] - thickness[0] / 2,
            (size[2] - top_length[0]) / 2
        ]
        self.top_mesh = Cuboid(thickness[0], size[0], top_length[0],
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

        self.semantic = 'Nozzle'


class Enveloping_Nozzle(ConceptTemplate):
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

        mesh_position = [
            0, 
            size[1] / 2,
            0
        ]
        self.mesh = Rectangular_Ring(size[1], size[0], size[2],
                                     size[0] - thickness[0] * 2, size[2] - thickness[0] * 2,
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

        self.semantic = 'Nozzle'


class Regular_Cover(ConceptTemplate):
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

        top_mesh_position = [
            0,
            (outer_size[1] + inner_size[1]) / 2,
            -outer_size[2] / 2
        ]
        self.top_mesh = Cuboid(outer_size[1] - inner_size[1], outer_size[0], outer_size[2],
                               position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0, 
            inner_size[1] / 2,
            -outer_size[2] / 2
        ]
        self.bottom_mesh = Rectangular_Ring(inner_size[1], outer_size[0], outer_size[2], inner_size[0], inner_size[2],
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

class Standard_Wick(ConceptTemplate):
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

        self.semantic = 'Wick'