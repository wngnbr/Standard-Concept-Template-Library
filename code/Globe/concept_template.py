import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh

class Standard_Sphere(ConceptTemplate):
    def __init__(self, radius, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.radius = radius

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.mesh = Sphere(radius[0])
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Sphere'


class Semi_Ring_Bracket(ConceptTemplate):
    def __init__(self, pivot_size, pivot_continuity, pivot_seperation, has_top_endpoint, has_bottom_endpoint, endpoint_radius, bracket_size, bracket_exist_angle, bracket_offset, bracket_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        bracket_exist_angle = [x / 180 * np.pi for x in bracket_exist_angle]
        bracket_rotation = [x / 180 * np.pi for x in bracket_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.pivot_size = pivot_size
        self.pivot_continuity = pivot_continuity
        self.pivot_seperation = pivot_seperation
        self.has_top_endpoint = has_top_endpoint
        self.has_bottom_endpoint = has_bottom_endpoint
        self.endpoint_radius = endpoint_radius
        self.bracket_size = bracket_size
        self.bracket_exist_angle = bracket_exist_angle
        self.bracket_offset = bracket_offset
        self.bracket_rotation = bracket_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        if pivot_continuity[0] == 1:
            self.pivot_mesh = Cylinder(pivot_size[1], pivot_size[0])
            vertices_list.append(self.pivot_mesh.vertices)
            faces_list.append(self.pivot_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.pivot_mesh.vertices)

        else:
            top_pivot_mesh_position = [
                0,
                pivot_seperation[0] / 2 + pivot_size[1] / 2,
                0
            ]
            self.top_pivot_mesh = Cylinder(pivot_size[1], pivot_size[0],
                                           position = top_pivot_mesh_position)
            vertices_list.append(self.top_pivot_mesh.vertices)
            faces_list.append(self.top_pivot_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.top_pivot_mesh.vertices)

            bottom_pivot_mesh_position = [
                0,
                -pivot_seperation[0] / 2 - pivot_size[1] / 2,
                0
            ]
            self.bottom_pivot_mesh = Cylinder(pivot_size[1], pivot_size[0],
                                              position = bottom_pivot_mesh_position)
            vertices_list.append(self.bottom_pivot_mesh.vertices)
            faces_list.append(self.bottom_pivot_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.bottom_pivot_mesh.vertices)

        if has_top_endpoint[0] == 1:
            if pivot_continuity[0] == 0:
                top_endpoint_mesh_position = [
                    0,
                    pivot_seperation[0] / 2 + pivot_size[1] + endpoint_radius[0],
                    0
                ]
            else:
                top_endpoint_mesh_position = [
                    0,
                    pivot_size[1] / 2 + endpoint_radius[0],
                    0
                ]
            self.top_endpoint_mesh = Sphere(endpoint_radius[0],
                                            position = top_endpoint_mesh_position)
            vertices_list.append(self.top_endpoint_mesh.vertices)
            faces_list.append(self.top_endpoint_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.top_endpoint_mesh.vertices)

        if has_bottom_endpoint[0] == 1:
            if pivot_continuity[0] == 0:
                bottom_endpoint_mesh_position = [
                    0,
                    -pivot_seperation[0] / 2 - pivot_size[1] - endpoint_radius[0],
                    0
                ]
            else:
                bottom_endpoint_mesh_position = [
                    0,
                    -pivot_size[1] / 2 - endpoint_radius[0],
                    0
                ]
            self.bottom_endpoint_mesh = Sphere(endpoint_radius[0],
                                               position = bottom_endpoint_mesh_position)
            vertices_list.append(self.bottom_endpoint_mesh.vertices)
            faces_list.append(self.bottom_endpoint_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.bottom_endpoint_mesh.vertices)

        bracket_mesh_position = [
            0,
            bracket_offset[0],
            0
        ]
        bracket_mesh_rotation = [
            0,
            -np.pi / 2 + bracket_exist_angle[0] / 2 - bracket_rotation[0],
            np.pi / 2
        ]
        self.bracket_mesh = Ring(bracket_size[2], bracket_size[0], bracket_size[1], bracket_exist_angle[0],
                                 position = bracket_mesh_position,
                                 rotation = bracket_mesh_rotation)
        vertices_list.append(self.bracket_mesh.vertices)
        faces_list.append(self.bracket_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bracket_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Bracket'


class Tilted_Bracket(ConceptTemplate):
    def __init__(self, pivot_size, bracket_size, circle_thickness, circle_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        circle_rotation = [x / 180 * np.pi for x in circle_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.pivot_size = pivot_size
        self.bracket_size = bracket_size
        self.circle_thickness = circle_thickness
        self.circle_rotation = circle_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        pivot_mesh_rotation = [np.pi / 2, 0, 0]
        self.pivot_mesh = Cylinder(pivot_size[1], pivot_size[0],
                                   rotation = pivot_mesh_rotation)
        vertices_list.append(self.pivot_mesh.vertices)
        faces_list.append(self.pivot_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.pivot_mesh.vertices)

        bracket_mesh_rotation = [
            np.pi / 2,
            np.pi / 2,
            0
        ]
        self.bracket_mesh = Ring(bracket_size[2], bracket_size[0], bracket_size[1], np.pi,
                                 rotation = bracket_mesh_rotation)
        vertices_list.append(self.bracket_mesh.vertices)
        faces_list.append(self.bracket_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bracket_mesh.vertices)

        circle_mesh_rotation = [
            0,
            0,
            circle_rotation[0]
        ]
        self.circle_mesh = Ring(circle_thickness[1], bracket_size[1], bracket_size[1] - circle_thickness[0],
                                rotation = circle_mesh_rotation)
        vertices_list.append(self.circle_mesh.vertices)
        faces_list.append(self.circle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.circle_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Bracket'


class Enclosed_Bracket(ConceptTemplate):
    def __init__(self, bracket_size, circle_radius, circle_thickness, half_circle_number, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.bracket_size = bracket_size
        self.circle_radius = circle_radius
        self.circle_thickness = circle_thickness
        self.half_circle_number = half_circle_number

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        bracket_mesh_rotation = [
            np.pi / 2,
            np.pi / 2,
            0
        ]
        self.bracket_mesh = Ring(bracket_size[2], bracket_size[0], bracket_size[1], np.pi, 
                                 rotation = bracket_mesh_rotation)
        vertices_list.append(self.bracket_mesh.vertices)
        faces_list.append(self.bracket_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.bracket_mesh.vertices)

        if half_circle_number[0] == 2:
            bracket_mesh_rotation = [
                np.pi / 2,
                0,
                0
            ]
            self.bracket_mesh = Ring(bracket_size[2], bracket_size[0], bracket_size[1], np.pi, 
                                     rotation = bracket_mesh_rotation)
            vertices_list.append(self.bracket_mesh.vertices)
            faces_list.append(self.bracket_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.bracket_mesh.vertices)

        outer_radius = circle_radius[0] + circle_thickness[0] / 2
        inner_radius = circle_radius[0] - circle_thickness[0] / 2
        self.circle_mesh = Ring(circle_thickness[1], outer_radius, inner_radius)
        vertices_list.append(self.circle_mesh.vertices)
        faces_list.append(self.circle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.circle_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Bracket'


class Cylindrical_Base(ConceptTemplate):
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

        top_mesh_position = [
            0,
            -top_size[2] / 2,
            0
        ]
        self.top_mesh = Cylinder(top_size[2], top_size[0], top_size[1],
                                 position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0,
            -top_size[2] - bottom_size[2] / 2,
            0
        ]
        self.bottom_mesh = Cylinder(bottom_size[2], bottom_size[0], bottom_size[1],
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

        self.semantic = 'Base'


class Cuboidal_Base(ConceptTemplate):
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

        top_mesh_position = [
            0,
            -top_size[1] / 2,
            0
        ]
        self.top_mesh = Cuboid(top_size[1], top_size[0], top_size[2],
                               position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0,
            -top_size[1] - bottom_size[1] / 2,
            0
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

        self.semantic = 'Base'


class Star_Shaped_Base(ConceptTemplate):
    def __init__(self, top_size, sub_size, num_legs, tilt_angle, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        tilt_angle = [x / 180 * np.pi for x in tilt_angle]
        super().__init__(position, rotation)

        # Record Parameters
        self.top_size = top_size
        self.sub_size = sub_size
        self.num_legs = num_legs
        self.tilt_angle = tilt_angle

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        top_mesh_position = [
            0,
            -top_size[1] / 2,
            0
        ]
        self.top_mesh = Cylinder(top_size[1], top_size[0], 
                                 position = top_mesh_position)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        for i in range(num_legs[0]):

            rot = np.pi * 2 / num_legs[0] * i
            claw_mesh_rotation = [
                tilt_angle[0],
                -rot,
                0
            ]
            x_ = sub_size[2] / 2 * np.cos(tilt_angle[0]) * np.sin(rot)
            z_ = sub_size[2] / 2 * np.cos(tilt_angle[0]) * np.cos(rot)
            y__ = -top_size[1] + sub_size[1] / 2 - sub_size[2] * np.sin(tilt_angle[0]) / 2
            claw_mesh_position = [
                -x_,
                y__,
                z_
            ]
            self.claw_mesh = Cuboid(sub_size[1], sub_size[0], sub_size[2],
                                    position = claw_mesh_position,
                                    rotation = claw_mesh_rotation)
            vertices_list.append(self.claw_mesh.vertices)
            faces_list.append(self.claw_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.claw_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Base'


class Special_Base(ConceptTemplate):
    def __init__(self, radius, top_size, top_rotation, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        top_rotation = [x / 180 * np.pi for x in top_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.radius = radius
        self.top_size = top_size
        self.top_rotation = top_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        top_mesh_position = [
            0,
            -top_size[1] / 2 * np.cos(top_rotation[0]),
            -top_size[1] / 2 * np.sin(top_rotation[0])
        ]
        top_mesh_rotation = [
            top_rotation[0],
            0,
            0
        ]
        self.top_mesh = Cylinder(top_size[1], top_size[0],
                                 position = top_mesh_position,
                                 rotation = top_mesh_rotation)
        vertices_list.append(self.top_mesh.vertices)
        faces_list.append(self.top_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.top_mesh.vertices)

        bottom_mesh_position = [
            0,
            -top_size[1] * np.cos(top_rotation[0]),
            radius[0] - top_size[1] * np.sin(top_rotation[0])
        ]
        self.bottom_mesh = Torus(radius[0], radius[1], 
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

        self.semantic = 'Base'


class Table_Like_Base(ConceptTemplate):
    def __init__(self, circle_size, num_legs, leg_size, leg_seperation, has_bottom_part, bottom_size, bottom_offset, position = [0, 0, 0], rotation = [0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.circle_size = circle_size
        self.num_legs = num_legs
        self.leg_size = leg_size
        self.leg_seperation = leg_seperation
        self.has_bottom_part = has_bottom_part
        self.bottom_size = bottom_size
        self.bottom_offset = bottom_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.ring_mesh = Ring(circle_size[2], circle_size[0], circle_size[1])
        vertices_list.append(self.ring_mesh.vertices)
        faces_list.append(self.ring_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.ring_mesh.vertices)

        for i in range(num_legs[0]):
            rotation_tmp = np.pi * 2 / num_legs[0] * i
            leg_mesh_position = [
                leg_seperation[0] * np.cos(rotation_tmp),
                -circle_size[2] / 2 - leg_size[1] / 2,
                leg_seperation[0] * np.sin(rotation_tmp)
            ]
            self.leg_mesh = Cylinder(leg_size[1], leg_size[0], 
                                     position = leg_mesh_position)
            vertices_list.append(self.leg_mesh.vertices)
            faces_list.append(self.leg_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.leg_mesh.vertices)

            if has_bottom_part[0] == 1:
                bottom_mesh_position = [
                    bottom_size[0] / 2 * np.cos(rotation_tmp),
                    -circle_size[2] / 2 - leg_size[1] + bottom_size[1] / 2 + bottom_offset[0],
                    bottom_size[0] / 2 * np.sin(rotation_tmp)
                ]
                bottom_mesh_rotation = [
                    0,
                    -rotation_tmp,
                    0
                ]
                self.bottom_mesh = Cuboid(bottom_size[1], bottom_size[0], bottom_size[2],
                                          position = bottom_mesh_position,
                                          rotation = bottom_mesh_rotation)
                vertices_list.append(self.bottom_mesh.vertices)
                faces_list.append(self.bottom_mesh.faces + total_num_vertices)
                total_num_vertices += len(self.bottom_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Base'