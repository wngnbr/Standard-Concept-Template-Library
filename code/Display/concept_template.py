import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh


class Frustum_Screen(ConceptTemplate):
    def __init__(self, has_additional_layer, additional_layer_size, size, back_part_offset, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.has_additional_layer = has_additional_layer
        self.additional_layer_size = additional_layer_size
        self.size = size
        self.back_part_offset = back_part_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        back_mesh_position = [0, 0, -size[2] / 2]
        back_mesh_rotation = [-np.pi / 2, 0, 0]
        self.back_mesh = Cuboid(size[2], size[0], size[1],
                                additional_layer_size[0], additional_layer_size[1],
                                [back_part_offset[0], back_part_offset[1]],
                                position=back_mesh_position,
                                rotation=back_mesh_rotation)
        vertices_list.append(self.back_mesh.vertices)
        faces_list.append(self.back_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.back_mesh.vertices)

        if has_additional_layer[0]:
            front_mesh_position = [0, 0, additional_layer_size[2] / 2]
            self.front_mesh = Cuboid(additional_layer_size[1], additional_layer_size[0], additional_layer_size[2],
                                     position=front_mesh_position)
            vertices_list.append(self.front_mesh.vertices)
            faces_list.append(self.front_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.front_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Screen'


class Standard_Screen(ConceptTemplate):
    def __init__(self, has_additional_layer, size, additional_layer_size, additional_layer_offset, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.has_additional_layer = has_additional_layer
        self.size = size
        self.additional_layer_size = additional_layer_size
        self.additional_layer_offset = additional_layer_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.front_mesh = Cuboid(size[1], size[0], size[2])
        vertices_list.append(self.front_mesh.vertices)
        faces_list.append(self.front_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.front_mesh.vertices)

        if has_additional_layer[0]:
            back_mesh_position = [
                additional_layer_offset[0],
                additional_layer_offset[1],
                -size[2] / 2 - additional_layer_size[2] / 2,
            ]
            self.back_mesh = Cuboid(additional_layer_size[1], additional_layer_size[0], additional_layer_size[2],
                                    position=back_mesh_position)
            vertices_list.append(self.back_mesh.vertices)
            faces_list.append(self.back_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.back_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Screen'


class Cuboidal_Base(ConceptTemplate):
    def __init__(self, size, base_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        base_rotation = [x / 180 * np.pi for x in base_rotation]
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.base_rotation = base_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        base_mesh_position = [
            0, 
            -size[1] / 2, 
            0
        ]
        base_mesh_rotation = [base_rotation[0], 0, 0]
        self.base_mesh = Cuboid(size[1], size[0], size[2],
                                position=base_mesh_position,
                                rotation=base_mesh_rotation)
        vertices_list.append(self.base_mesh.vertices)
        faces_list.append(self.base_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.base_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Base'


class Round_Base(ConceptTemplate):
    def __init__(self, size, base_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        base_rotation = [x / 180 * np.pi for x in base_rotation]
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.base_rotation = base_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # meshes definition
        base_mesh_position = [
            0, 
            -size[1] / 2, 
            0
        ]
        base_mesh_rotation = [base_rotation[0], 0, 0]
        self.base_mesh = Cylinder(size[1], size[0], size[0],
                                  position=base_mesh_position,
                                  rotation=base_mesh_rotation)
        vertices_list.append(self.base_mesh.vertices)
        faces_list.append(self.base_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.base_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Base'


class TShaped_Base(ConceptTemplate):
    def __init__(self, main_size, sub_size, base_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        base_rotation = [x / 180 * np.pi for x in base_rotation]
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.main_size = main_size
        self.sub_size = sub_size
        self.base_rotation = base_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        sub_mesh_position = [
            0, 
            -main_size[1] / 2, 
            0
        ]
        sub_mesh_rotation = [base_rotation[0], 0, 0]
        self.sub_mesh = Cuboid(main_size[1], sub_size[0], sub_size[1],
                               position=sub_mesh_position,
                               rotation=sub_mesh_rotation)
        vertices_list.append(self.sub_mesh.vertices)
        faces_list.append(self.sub_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.sub_mesh.vertices)

        main_mesh_position = [
            0,
            -main_size[1] / 2 - (main_size[2] + sub_size[1]) * np.sin(base_rotation[0]) / 2,
            (main_size[2] + sub_size[1]) * np.cos(base_rotation[0]) / 2,
        ]
        main_mesh_rotation = [base_rotation[0], 0, 0]
        self.main_mesh = Cuboid(main_size[1], main_size[0], main_size[2],
                                position=main_mesh_position,
                                rotation=main_mesh_rotation)
        vertices_list.append(self.main_mesh.vertices)
        faces_list.append(self.main_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.main_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Base'


class Vshaped_Base(ConceptTemplate):
    def __init__(self, size, open_angle, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        open_angle = [x / 180 * np.pi for x in open_angle]
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.open_angle = open_angle

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        left_mesh_position = [0, -size[1] / 2, size[2] / 2]
        self.left_mesh = Cuboid(size[1], size[0], size[2], 
                                position=left_mesh_position)
        self.left_mesh.vertices = apply_transformation(self.left_mesh.vertices, position=[0, 0, 0], rotation=[0, -open_angle[0] / 2, 0])
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [0, -size[1] / 2, size[2] / 2]
        self.right_mesh = Cuboid(size[1], size[0], size[2], 
                                 position=right_mesh_position)
        self.right_mesh.vertices = apply_transformation(self.right_mesh.vertices, position=[0, 0, 0], rotation=[0, open_angle[0] / 2, 0])
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Base'


class CuboidalRear_Support(ConceptTemplate):
    def __init__(self, number_of_supports, size, separation, support_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        support_rotation = [x / 180 * np.pi for x in support_rotation]
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_supports = number_of_supports
        self.size = size
        self.separation = separation
        self.support_rotation = support_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(number_of_supports[0]):
            mesh_position = [
                (separation[0] + size[0]) * i, 
                -size[1] / 2, 
                -size[2] / 2
            ]
            self.mesh = Cuboid(size[1], size[0], size[2], 
                               position=mesh_position)
            self.mesh.vertices = apply_transformation(self.mesh.vertices, position=[0, 0, 0], rotation=[support_rotation[0], 0, 0])
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Support'


class Cuboidal_Support(ConceptTemplate):
    def __init__(self, number_of_supports, size, separation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_supports = number_of_supports
        self.size = size
        self.separation = separation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(number_of_supports[0]):
            mesh_position = [
                (separation[0] + size[0]) * i,
                -size[1] / 2,
                0,
            ]
            self.mesh = Cuboid(size[1], size[0], size[2], 
                               position=mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Support'


class Trifold_Support(ConceptTemplate):
    def __init__(self, has_upper_part, upper_part_size, upper_offset, has_middle_part, middle_part_size, middle_offset, has_bottom_part, bottom_part_size, rotation_of_bottom, bottom_offset, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation_of_bottom = [x / 180 * np.pi for x in rotation_of_bottom]
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.has_upper_part = has_upper_part
        self.upper_part_size = upper_part_size
        self.upper_offset = upper_offset
        self.has_middle_part = has_middle_part
        self.middle_part_size = middle_part_size
        self.middle_offset = middle_offset
        self.has_bottom_part = has_bottom_part
        self.bottom_part_size = bottom_part_size
        self.rotation_of_bottom = rotation_of_bottom
        self.bottom_offset = bottom_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        base_position = [0, 0, 0]

        if has_upper_part[0]:
            mesh_tmp_position = [
                upper_offset[0],
                upper_offset[1],
                upper_offset[2] - upper_part_size[2] / 2
            ]
            self.mesh_tmp = Cuboid(upper_part_size[1], upper_part_size[0], upper_part_size[2],
                                  position=mesh_tmp_position)
            vertices_list.append(self.mesh_tmp.vertices)
            faces_list.append(self.mesh_tmp.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_tmp.vertices)
            base_position = [
                base_position[0] + upper_offset[0],
                base_position[1] + upper_offset[1],
                base_position[2] + upper_offset[2] - upper_part_size[2],
            ]

        if has_middle_part[0]:
            mesh_tmp_position = [
                base_position[0],
                base_position[1] + middle_offset[0] - middle_part_size[1] / 2,
                base_position[2] - middle_part_size[2] / 2
            ]
            self.mesh_tmp = Cuboid(middle_part_size[1], middle_part_size[0], middle_part_size[2],
                                   position=mesh_tmp_position)
            vertices_list.append(self.mesh_tmp.vertices)
            faces_list.append(self.mesh_tmp.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_tmp.vertices)
            base_position = [
                base_position[0],
                base_position[1] + middle_offset[0] - middle_part_size[1],
                base_position[2] - middle_part_size[2],
            ]

        if has_bottom_part[0]:
            base_position = [
                base_position[0],
                base_position[1],
                base_position[2] + bottom_offset[0],
            ]
            mesh_tmp_position = [
                0, 
                -bottom_part_size[1] / 2, 
                bottom_part_size[2] / 2
            ]
            self.mesh_tmp = Cuboid(bottom_part_size[1], bottom_part_size[0], bottom_part_size[2],
                                   position=mesh_tmp_position)
            self.mesh_tmp.vertices = apply_transformation(self.mesh_tmp.vertices, base_position, [rotation_of_bottom[0], 0, 0])
            vertices_list.append(self.mesh_tmp.vertices)
            faces_list.append(self.mesh_tmp.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_tmp.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Support'


class TShaped_Support(ConceptTemplate):
    def __init__(self, upper_part_size, upper_offset, middle_part_size, middle_offset, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.upper_part_size = upper_part_size
        self.upper_offset = upper_offset
        self.middle_part_size = middle_part_size
        self.middle_offset = middle_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        top_mesh_position = [
            upper_offset[0], 
            upper_offset[1], 
            upper_offset[2] - upper_part_size[2] / 2
        ]
        self.top_mesh = Cuboid(upper_part_size[1], upper_part_size[0], upper_part_size[2],
                               position=top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        middle_mesh_position = [
            upper_offset[0],
            upper_offset[1] - middle_part_size[1] / 2 - upper_part_size[1] / 2,
            upper_offset[2] - middle_part_size[2] / 2 + middle_offset[0],
        ]
        self.middle_mesh = Cuboid(middle_part_size[1], middle_part_size[0], middle_part_size[2],
                                  position=middle_mesh_position)
        vertices_list.append(self.middle_mesh.vertices)
        faces_list.append(self.middle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.middle_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Support'