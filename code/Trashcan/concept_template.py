import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation, adjust_position_from_rotation, list_add
from knowledge_utils import *
import trimesh

class Cylindrical_Body(ConceptTemplate):
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

        middle_radius = outer_size[0] * (1 - inner_size[2] / outer_size[2]) + outer_size[1] * inner_size[2] / outer_size[2]
        bottom_height = outer_size[2] - inner_size[2]
        bottom_mesh_position = [0, -inner_size[2] / 2, 0]
        self.bottom_mesh = Cylinder(bottom_height, middle_radius, outer_size[1], 
                                    position=bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [0, (outer_size[2] - inner_size[2]) / 2, 0]
        self.top_mesh = Ring(inner_size[2], outer_size[0], inner_size[0], 
                             outer_bottom_radius = middle_radius,
                             inner_bottom_radius = inner_size[1],
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

        self.semantic = 'Body'


class Prismatic_Body(ConceptTemplate):
    def __init__(self, top_size, bottom_size, height, top_offset, thickness, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.top_size = top_size
        self.bottom_size = bottom_size
        self.height = height
        self.top_offset = top_offset
        self.thickness = thickness

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0, 
            -(height[1] - height[0]) / 2, 
            0
        ]
        self.mesh = Rectangular_Ring(height[0], top_size[0], top_size[1], 
                                    top_size[0] - thickness[0] * 2, top_size[1] - thickness[1] * 2, 
                                    outer_bottom_length = bottom_size[0], outer_bottom_width = bottom_size[1],
                                    inner_bottom_length = bottom_size[0] - thickness[0] * 2, inner_bottom_width = bottom_size[1] - thickness[1] * 2,
                                    back_height = height[1],
                                    top_bottom_offset = [top_offset[0], top_offset[1]],
                                    position = mesh_position)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        middle_x = (top_size[0] - thickness[0] * 2) * thickness[2] / height[0] + (bottom_size[0] - thickness[0] * 2) * (height[0] - thickness[2]) / height[0]
        middle_z = (top_size[1] - thickness[1] * 2) * thickness[2] / height[0] + (bottom_size[1] - thickness[1] * 2) * (height[0] - thickness[2]) / height[0]
        middle_offset_x = top_offset[0] * thickness[2] / height[0]
        middle_offset_z = top_offset[1] * thickness[2] / height[0]
        bottom_mesh_position = [
            0, 
            -(height[1] - thickness[2]) / 2, 
            0
        ]
        self.bottom_mesh = Cuboid(thickness[2], middle_x, middle_z, 
                                  bottom_length = bottom_size[0] - thickness[0] * 2, bottom_width = bottom_size[1] - thickness[1] * 2,
                                  top_offset = [middle_offset_x, middle_offset_z],
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

        self.semantic = 'Body'

        
class Separated_Cylindrical_Body(ConceptTemplate):
    def __init__(self, outer_size, inner_size, clapboard_size, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.outer_size = outer_size
        self.inner_size = inner_size
        self.clapboard_size = clapboard_size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        middle_radius = outer_size[0] * (1 - inner_size[2] / outer_size[2]) + outer_size[1] * inner_size[2] / outer_size[2]
        bottom_height = outer_size[2] - inner_size[2]
        bottom_mesh_position = [0, -inner_size[2] / 2, 0]
        self.bottom_mesh = Cylinder(bottom_height, middle_radius, outer_size[1], 
                                    position=bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [0, (outer_size[2] - inner_size[2]) / 2, 0]
        self.top_mesh = Ring(inner_size[2], outer_size[0], inner_size[0], 
                             outer_bottom_radius = middle_radius,
                             inner_bottom_radius = inner_size[1],
                             position=top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        clapboard_mesh_position = [
            0, 
            clapboard_size[1] / 2 - (outer_size[2] / 2 - (outer_size[2] - inner_size[2])), 
            0
        ]
        self.clapboard_mesh = Cuboid(clapboard_size[1], clapboard_size[0], clapboard_size[2], 
                                     position = clapboard_mesh_position)
        vertices_list.append(self.clapboard_mesh.vertices)
        faces_list.append(self.clapboard_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.clapboard_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Cylindrical_Cover(ConceptTemplate):
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

        self.semantic = 'Cover'


class Cuboidal_Cover(ConceptTemplate):
    def __init__(self, top_size, bottom_size, height, top_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.top_size = top_size
        self.bottom_size = bottom_size
        self.height = height
        self.top_offset = top_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0,
            height[0] / 2,
            bottom_size[1] / 2
        ]
        self.mesh = Cuboid(height[0], top_size[0], top_size[1], 
                           bottom_size[0], bottom_size[1], 
                           [top_offset[0], top_offset[1]],
                           back_height = height[1],
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

        self.semantic = 'Cover'


class Double_Layer_Cuboidal_Cover(ConceptTemplate):
    def __init__(self, top_size, bottom_size, top_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.top_size = top_size
        self.bottom_size = bottom_size
        self.top_offset = top_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bottom_mesh_position = [
            0,
            bottom_size[1] / 2,
            bottom_size[2] / 2
        ]
        self.bottom_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                                  position = bottom_mesh_position)
        vertices_list.append(self.bottom_mesh.vertices)
        faces_list.append(self.bottom_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bottom_mesh.vertices)

        top_mesh_position = [
            top_offset[0],
            bottom_size[1] + top_size[1] / 2,
            bottom_size[2] / 2 + top_offset[1]
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

        self.semantic = 'Cover'


class Cylindrical_Hollow_Cover(ConceptTemplate):
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
            size[2] / 2,
            0
        ]
        self.mesh = Ring(size[2], size[0], size[1],
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

        self.semantic = 'Cover'


class Cuboidal_Hollow_Cover(ConceptTemplate):
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

        mesh_position = [
            0,
            outer_size[1] / 2,
            0
        ]
        self.mesh = Rectangular_Ring(outer_size[1], outer_size[0], outer_size[2],
                                    inner_size[0], inner_size[1], 
                                    [inner_offset[0], inner_offset[1]],
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

        self.semantic = 'Cover'


class Holed_Cylindrical_Cover(ConceptTemplate):
    def __init__(self, radius, height, exist_angle, num_sides, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        exist_angle = [x / 180 * np.pi for x in exist_angle]
        super().__init__(position, rotation)

        # Record Parameters
        self.radius = radius
        self.height = height
        self.exist_angle = exist_angle
        self.num_sides = num_sides

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(num_sides[0]):
            rotation_y = i * np.pi * 2 / num_sides[0]
            mesh_position = [
                0,
                height[2] / 2,
                0
            ]
            mesh_rotation = [0, rotation_y, 0]
            self.mesh = Ring(height[2], radius[0], radius[1], exist_angle[0],
                             position = mesh_position,
                             rotation = mesh_rotation)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        top_mesh_position = [
            0,
            height[2] + height[1] + (height[0] - height[1]) / 2,
            0
        ]
        self.top_mesh = Cylinder(height[0] - height[1], radius[0], 
                                 position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        middle_mesh_position = [
            0,
            height[2] + height[1] / 2,
            0
        ]
        self.middle_mesh = Ring(height[1], radius[0], radius[1], 
                                position = middle_mesh_position)
        vertices_list.append(self.middle_mesh.vertices)
        faces_list.append(self.middle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.middle_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Cover'


class Holed_Cuboidal_Cover(ConceptTemplate):
    def __init__(self, outer_size, inner_size, front_behind_hole_size, left_right_hole_size, has_hole, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.outer_size = outer_size
        self.inner_size = inner_size
        self.front_behind_hole_size = front_behind_hole_size
        self.left_right_hole_size = left_right_hole_size
        self.has_hole = has_hole

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        # top part
        top_mesh_position = [
            0,
            outer_size[1] / 2 + inner_size[1] / 2,
            0
        ]
        self.top_mesh = Cuboid(outer_size[1] - inner_size[1], outer_size[0], outer_size[2], 
                               position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        # front part
        if has_hole[0] == 0:
            mesh_position = [
                0,
                inner_size[1] / 2,
                (outer_size[2] + inner_size[2]) / 4
            ]
            self.mesh = Cuboid(inner_size[1], outer_size[0], (outer_size[2] - inner_size[2]) / 2, 
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        elif has_hole[0] == 1:
            mesh_1_position = [
                -(outer_size[0] + front_behind_hole_size[0]) / 4,
                inner_size[1] / 2,
                (outer_size[2] + inner_size[2]) / 4
            ]
            self.mesh_1 = Cuboid(inner_size[1], (outer_size[0] - front_behind_hole_size[0]) / 2, (outer_size[2] - inner_size[2]) / 2, 
                                 position = mesh_1_position)
            vertices_list.append(self.mesh_1.vertices)
            faces_list.append(self.mesh_1.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_1.vertices)

            mesh_2_position = [
                (outer_size[0] + front_behind_hole_size[0]) / 4,
                inner_size[1] / 2,
                (outer_size[2] + inner_size[2]) / 4
            ]
            self.mesh_2 = Cuboid(inner_size[1], (outer_size[0] - front_behind_hole_size[0]) / 2, (outer_size[2] - inner_size[2]) / 2, 
                                 position = mesh_2_position)
            vertices_list.append(self.mesh_2.vertices)
            faces_list.append(self.mesh_2.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_2.vertices)

            mesh_3_position = [
                0,
                (inner_size[1] + front_behind_hole_size[1]) / 2,
                (outer_size[2] + inner_size[2]) / 4
            ]
            self.mesh_3 = Cuboid(inner_size[1] - front_behind_hole_size[1], front_behind_hole_size[0], (outer_size[2] - inner_size[2]) / 2, 
                                 position = mesh_3_position)
            vertices_list.append(self.mesh_3.vertices)
            faces_list.append(self.mesh_3.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_3.vertices)

        # behind part
        if has_hole[1] == 0:
            mesh_position = [
                0,
                inner_size[1] / 2,
                -(outer_size[2] + inner_size[2]) / 4
            ]
            self.mesh = Cuboid(inner_size[1], outer_size[0], (outer_size[2] - inner_size[2]) / 2, 
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        elif has_hole[1] == 1:
            mesh_1_position = [
                -(outer_size[0] + front_behind_hole_size[0]) / 4,
                inner_size[1] / 2,
                -(outer_size[2] + inner_size[2]) / 4
            ]
            self.mesh_1 = Cuboid(inner_size[1], (outer_size[0] - front_behind_hole_size[0]) / 2, (outer_size[2] - inner_size[2]) / 2, 
                                 position = mesh_1_position)
            vertices_list.append(self.mesh_1.vertices)
            faces_list.append(self.mesh_1.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_1.vertices)

            mesh_2_position = [
                (outer_size[0] + front_behind_hole_size[0]) / 4,
                inner_size[1] / 2,
                -(outer_size[2] + inner_size[2]) / 4
            ]
            self.mesh_2 = Cuboid(inner_size[1], (outer_size[0] - front_behind_hole_size[0]) / 2, (outer_size[2] - inner_size[2]) / 2, 
                                 position = mesh_2_position)
            vertices_list.append(self.mesh_2.vertices)
            faces_list.append(self.mesh_2.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_2.vertices)

            mesh_3_position = [
                0,
                (inner_size[1] + front_behind_hole_size[1]) / 2,
                -(outer_size[2] + inner_size[2]) / 4
            ]
            self.mesh_3 = Cuboid(inner_size[1] - front_behind_hole_size[1], front_behind_hole_size[0], (outer_size[2] - inner_size[2]) / 2, 
                                 position = mesh_3_position)
            vertices_list.append(self.mesh_3.vertices)
            faces_list.append(self.mesh_3.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_3.vertices)

        # left part
        if has_hole[2] == 0:
            mesh_position = [
                -(outer_size[0] + inner_size[0]) / 4,
                inner_size[1] / 2,
                0
            ]
            self.mesh = Cuboid(inner_size[1], (outer_size[0] - inner_size[0]) / 2, inner_size[2], 
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        elif has_hole[2] == 1:
            mesh_1_position = [
                -(outer_size[0] + inner_size[0]) / 4,
                inner_size[1] / 2,
                -(inner_size[2] + left_right_hole_size[0]) / 4
            ]
            self.mesh_1 = Cuboid(inner_size[1], (outer_size[0] - inner_size[0]) / 2, (inner_size[2] - left_right_hole_size[0]) / 2,  
                                 position = mesh_1_position)
            vertices_list.append(self.mesh_1.vertices)
            faces_list.append(self.mesh_1.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_1.vertices)

            mesh_2_position = [
                -(outer_size[0] + inner_size[0]) / 4,
                inner_size[1] / 2,
                (inner_size[2] + left_right_hole_size[0]) / 4
            ]
            self.mesh_2 = Cuboid(inner_size[1], (outer_size[0] - inner_size[0]) / 2, (inner_size[2] - left_right_hole_size[0]) / 2,  
                                 position = mesh_2_position)
            vertices_list.append(self.mesh_2.vertices)
            faces_list.append(self.mesh_2.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_2.vertices)

            mesh_3_position = [
                -(outer_size[0] + inner_size[0]) / 4,
                (inner_size[1] + left_right_hole_size[1]) / 2,
                0
            ]
            self.mesh_3 = Cuboid(inner_size[1] - left_right_hole_size[1], (outer_size[2] - inner_size[2]) / 2, left_right_hole_size[0], 
                                 position = mesh_3_position)
            vertices_list.append(self.mesh_3.vertices)
            faces_list.append(self.mesh_3.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_3.vertices)

        # right part
        if has_hole[3] == 0:
            mesh_position = [
                (outer_size[0] + inner_size[0]) / 4,
                inner_size[1] / 2,
                0
            ]
            self.mesh = Cuboid(inner_size[1], (outer_size[0] - inner_size[0]) / 2, inner_size[2], 
                               position = mesh_position)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        elif has_hole[3] == 1:
            mesh_1_position = [
                (outer_size[0] + inner_size[0]) / 4,
                inner_size[1] / 2,
                -(inner_size[2] + left_right_hole_size[0]) / 4
            ]
            self.mesh_1 = Cuboid(inner_size[1], (outer_size[0] - inner_size[0]) / 2, (inner_size[2] - left_right_hole_size[0]) / 2,  
                                 position = mesh_1_position)
            vertices_list.append(self.mesh_1.vertices)
            faces_list.append(self.mesh_1.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_1.vertices)

            mesh_2_position = [
                (outer_size[0] + inner_size[0]) / 4,
                inner_size[1] / 2,
                (inner_size[2] + left_right_hole_size[0]) / 4
            ]
            self.mesh_2 = Cuboid(inner_size[1], (outer_size[0] - inner_size[0]) / 2, (inner_size[2] - left_right_hole_size[0]) / 2,  
                                 position = mesh_2_position)
            vertices_list.append(self.mesh_2.vertices)
            faces_list.append(self.mesh_2.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_2.vertices)

            mesh_3_position = [
                (outer_size[0] + inner_size[0]) / 4,
                (inner_size[1] + left_right_hole_size[1]) / 2,
                0
            ]
            self.mesh_3 = Cuboid(inner_size[1] - left_right_hole_size[1], (outer_size[2] - inner_size[2]) / 2, left_right_hole_size[0], 
                                 position = mesh_3_position)
            vertices_list.append(self.mesh_3.vertices)
            faces_list.append(self.mesh_3.faces + total_num_vertices)
            total_num_vertices += len(self.mesh_3.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Cover'


class Standard_Wheel(ConceptTemplate):
    def __init__(self, size, seperation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.seperation = seperation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        left_mesh_position = [
            -seperation[0],
            0,
            0
        ]
        left_mesh_rotation = [0, 0, np.pi / 2]
        self.left_mesh = Cylinder(size[1], size[0], 
                                  position = left_mesh_position,
                                  rotation = left_mesh_rotation)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            seperation[0],
            0,
            0
        ]
        right_mesh_rotation = [0, 0, np.pi / 2]
        self.right_mesh = Cylinder(size[1], size[0], 
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


class Cylindrical_Shell(ConceptTemplate):
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

        self.mesh = Ring(outer_size[2], outer_size[0], inner_size[0],
                         outer_bottom_radius = outer_size[1], inner_bottom_radius = inner_size[1])
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Shell'


class Cuboidal_Shell(ConceptTemplate):
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

        self.mesh = Rectangular_Ring(outer_size[1], outer_size[0], outer_size[2],
                                     inner_size[0], inner_size[1])
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Shell'