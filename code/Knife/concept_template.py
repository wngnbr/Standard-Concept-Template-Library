import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation, adjust_position_from_rotation, list_add
from knowledge_utils import *
import trimesh

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

        self.semantic = 'Handle'


class T_Shaped_Handle(ConceptTemplate):
    def __init__(self, main_size, bottom_size, bottom_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.main_size = main_size
        self.bottom_size = bottom_size
        self.bottom_offset = bottom_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.mesh = Cuboid(main_size[1], main_size[0], main_size[2])
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        bottom_mesh_position = [
            bottom_offset[0], 
            -(main_size[1] + bottom_size[1]) / 2,
            bottom_offset[1]
        ]
        self.bottom_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
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

        self.mesh = Cylinder(size[2], size[0], size[1])
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


class Curved_Handle(ConceptTemplate):
    def __init__(self, radius, thickness, exist_angle, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        exist_angle = [x / 180 * np.pi for x in exist_angle]
        super().__init__(position, rotation)

        # Record Parameters
        self.radius = radius
        self.thickness = thickness
        self.exist_angle = exist_angle

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0,
            0,
            (radius[0] + radius[1]) / 2
        ]
        mesh_rotation = [np.pi / 2, np.pi / 2, 0]
        self.mesh = Ring(thickness[0], radius[0], radius[1], exist_angle[0],
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


class Multideck_Handle(ConceptTemplate):
    def __init__(self, bottom_size, beside_size, beside_seperation, beside_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.bottom_size = bottom_size
        self.beside_size = beside_size
        self.beside_seperation = beside_seperation
        self.beside_offset = beside_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            0,
            -beside_size[1] / 2,
            0
        ]
        self.bottom_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        left_mesh_position = [
            beside_seperation[0] / 2 + beside_offset[0],
            bottom_size[1] / 2,
            beside_offset[1]
        ]
        self.left_mesh = Cuboid(beside_size[1], beside_size[0], beside_size[2],
                                position = left_mesh_position)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            -beside_seperation[0] / 2 + beside_offset[0],
            bottom_size[1] / 2,
            beside_offset[1]
        ]
        self.right_mesh = Cuboid(beside_size[1], beside_size[0], beside_size[2],
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

        self.semantic = 'Handle'


class Enveloping_Handle(ConceptTemplate):
    def __init__(self, size, thickness, gap_width, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.thickness = thickness
        self.gap_width = gap_width

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            0,
            -(size[1] - thickness[0]) / 2,
            0
        ]
        self.bottom_mesh = Cuboid(thickness[0], size[0], size[2],
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        behind_mesh_position = [
            -(size[0] - thickness[1]) / 2,
            thickness[0] / 2,
            0
        ]
        self.behind_mesh = Cuboid(size[1] - thickness[0], thickness[1], size[2],
                                  position = behind_mesh_position)
        vertices_list.append(self.behind_mesh.vertices)
        faces_list.append(self.behind_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.behind_mesh.vertices)

        left_mesh_position = [
            thickness[1] / 2,
            thickness[0] / 2,
            (size[2] - thickness[2]) / 2
        ]
        self.left_mesh = Cuboid(size[1] - thickness[0], size[0] - thickness[1], thickness[2],
                                position = left_mesh_position)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            thickness[1] / 2,
            thickness[0] / 2,
            -(size[2] - thickness[2]) / 2
        ]
        self.right_mesh = Cuboid(size[1] - thickness[0], size[0] - thickness[1], thickness[2],
                                 position = right_mesh_position)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        left_front_mesh_position = [
            (size[0] - thickness[2]) / 2,
            thickness[0] / 2,
            size[2] / 2 - thickness[2] / 2 - (size[2] - gap_width[0]) / 4
        ]
        self.left_front_mesh = Cuboid(size[1] - thickness[0], thickness[2], (size[2] - gap_width[0]) / 2 - thickness[2],
                                      position = left_front_mesh_position)
        vertices_list.append(self.left_front_mesh.vertices)
        faces_list.append(self.left_front_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_front_mesh.vertices)

        right_front_mesh_position = [
            (size[0] - thickness[2]) / 2,
            thickness[0] / 2,
            -size[2] / 2 + thickness[2] / 2 + (size[2] - gap_width[0]) / 4
        ]
        self.right_front_mesh = Cuboid(size[1] - thickness[0], thickness[2], (size[2] - gap_width[0]) / 2 - thickness[2],
                                       position = right_front_mesh_position)
        vertices_list.append(self.right_front_mesh.vertices)
        faces_list.append(self.right_front_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_front_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Handle'


class Cuboidal_Blade(ConceptTemplate):
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

        self.semantic = 'Blade'


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

        mesh_position = [
            0, 
            root_size[1] / 2,
            0
        ]
        self.mesh = Cuboid(root_size[1], root_size[0], root_size[3],
                           root_size[0], root_size[2],
                           top_offset = [0, root_z_offset[0]],
                           position = mesh_position)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        tip_mesh_position = [
            0, 
            root_size[1] + tip_length[0] / 2,
            root_z_offset[0]
        ]
        self.tip_mesh = Cuboid(tip_length[0], root_size[0], 0,
                               root_size[0], root_size[3],
                               top_offset = [0, tip_z_offset[0]],
                               position = tip_mesh_position)
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
    def __init__(self, root_size, root_z_offset, tip_length, tip_angle, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        tip_angle = [x / 180 * np.pi for x in tip_angle]
        super().__init__(position, rotation)

        # Record Parameters
        self.root_size = root_size
        self.root_z_offset = root_z_offset
        self.tip_length = tip_length
        self.tip_angle = tip_angle

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0, 
            root_size[1] / 2,
            0
        ]
        self.mesh = Cuboid(root_size[1], root_size[0], root_size[3],
                           root_size[0], root_size[2],
                           top_offset = [0, root_z_offset[0]],
                           position = mesh_position)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        tip_mesh_position = [
            0, 
            root_size[1],
            root_z_offset[0] - root_size[3] / 2
        ]
        tip_mesh_rotation = [0, 0, np.pi / 2]
        self.tip_mesh = Cylinder(root_size[0], tip_length[0], tip_length[0],
                                 root_size[3], root_size[3],
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


class Standard_Guard(ConceptTemplate):
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

        self.semantic = 'Guard'


class Regular_Button(ConceptTemplate):
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

        self.semantic = 'Button'