import copy
import open3d as o3d
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh


class Cylindrical_body(ConceptTemplate):
    def __init__(self, num_of_part, all_sizes, x_z_ratio, position=[0, 0, 0], rotation=[0, 0, 0]):
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.num_of_part = num_of_part
        self.all_sizes = [all_sizes[0:3]] + [[all_sizes[i], all_sizes[i + 1]] for i in range(3, len(all_sizes), 2)]
        self.x_z_ratio = x_z_ratio

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        top_radius, bottom_radius, height = self.all_sizes[0][0], self.all_sizes[0][2], self.all_sizes[0][1]
        self.mesh_1 = Cylinder(height=height, top_radius=top_radius, bottom_radius=bottom_radius, top_radius_z=top_radius * x_z_ratio[0],
                               bottom_radius_z=bottom_radius * x_z_ratio[0])
        vertices_list.append(self.mesh_1.vertices)
        faces_list.append(self.mesh_1.faces + total_num_vertices)
        total_num_vertices += len(self.mesh_1.vertices)

        delta_height = self.all_sizes[0][1] / 2
        for part_idx in range(1, num_of_part[0]):
            top_radius, bottom_radius, height = self.all_sizes[part_idx][0], self.all_sizes[part_idx - 1][0], self.all_sizes[part_idx][1]
            delta_height += height
            part_mesh_position = [0, delta_height - height / 2, 0]
            self.part_mesh = Cylinder(height=height, top_radius=top_radius, bottom_radius=bottom_radius, top_radius_z=top_radius * x_z_ratio[0],
                                      bottom_radius_z=bottom_radius * x_z_ratio[0],
                                      position=part_mesh_position)
            vertices_list.append(self.part_mesh.vertices)
            faces_list.append(self.part_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.part_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Cuboidal_body(ConceptTemplate):
    def __init__(self, size, position=[0, 0, 0], rotation=[0, 0, 0]):
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.mesh_1 = Cuboid(size[1], size[0], size[2])
        vertices_list.append(self.mesh_1.vertices)
        faces_list.append(self.mesh_1.faces + total_num_vertices)
        total_num_vertices += len(self.mesh_1.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Toothpaste_body(ConceptTemplate):
    def __init__(self, radius, bottom_length, height, position=[0, 0, 0], rotation=[0, 0, 0]):
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.radius = radius
        self.bottom_length = bottom_length
        self.height = height

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.base_mesh = Cylinder(height[0], radius[0],
                                  bottom_radius=1e-2,
                                  bottom_radius_z=bottom_length[0] / 2)
        vertices_list.append(self.base_mesh.vertices)
        faces_list.append(self.base_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.base_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Regular_nozzle(ConceptTemplate):
    def __init__(self, num_of_part, num_of_nozzle, nozzle_size, nozzle_length, nozzle_offset, nozzle_rotation, parts_params,
                 position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        nozzle_rotation = [x / 180 * np.pi for x in nozzle_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.num_of_part = num_of_part
        self.num_of_nozzle = num_of_nozzle
        self.nozzle_size = nozzle_size
        self.nozzle_length = nozzle_length
        self.nozzle_offset = nozzle_offset
        self.nozzle_rotation = nozzle_rotation
        self.parts_params = [[parts_params[i], parts_params[i + 1]] for i in range(0, len(parts_params), 2)]

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        delta_height = 0
        for part_idx in range(num_of_part[0]):
            delta_height += self.parts_params[part_idx][1]
            part_mesh_position = [0, delta_height - self.parts_params[part_idx][1] / 2, 0]
            self.part_mesh = Cylinder(self.parts_params[part_idx][1], self.parts_params[part_idx][0], position=part_mesh_position)
            vertices_list.append(self.part_mesh.vertices)
            faces_list.append(self.part_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.part_mesh.vertices)

        nozzle_mesh_position = [0, delta_height + nozzle_offset[0] - nozzle_length[0] / 2 * np.sin(nozzle_rotation[0]),
                                nozzle_length[0] / 2 * np.cos(nozzle_rotation[0]) + self.parts_params[num_of_part[0] - 1][0]]
        self.nozzle_mesh = Cuboid(nozzle_size[1], nozzle_size[0], nozzle_length[0], position=nozzle_mesh_position)
        self.nozzle_mesh.vertices = apply_transformation(self.nozzle_mesh.vertices, [0, 0, 0], [0, nozzle_rotation[2], 0])
        vertices_list.append(self.nozzle_mesh.vertices)
        faces_list.append(self.nozzle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.nozzle_mesh.vertices)

        if num_of_nozzle == 2:
            other_nozzle_mesh_position = [0,
                                          delta_height + nozzle_offset[0] - nozzle_length[0] / 2 * np.sin(nozzle_rotation[0]) - nozzle_length[1] / 2 * np.sin(nozzle_rotation[1]),
                                          nozzle_length[0] * np.cos(nozzle_rotation[0]) + self.parts_params[num_of_part[0] - 1][0] + nozzle_length[1] * np.cos(
                                              nozzle_rotation[1]) / 2]
            self.other_nozzle_mesh = Cuboid(nozzle_size[1], nozzle_size[0], nozzle_length[0], position=other_nozzle_mesh_position)
            self.other_nozzle_mesh.vertices = apply_transformation(self.other_nozzle_mesh.vertices, [0, 0, 0], [0, nozzle_rotation[2], 0])
            vertices_list.append(self.other_nozzle_mesh.vertices)
            faces_list.append(self.other_nozzle_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.other_nozzle_mesh.vertices)
        else:
            pass

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Nozzle'


class Cylindrical_cap(ConceptTemplate):
    def __init__(self, outer_size, inner_size, x_z_ratio, position=[0, 0, 0], rotation=[0, 0, 0]):
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.outer_size = outer_size
        self.inner_size = inner_size
        self.x_z_ratio = x_z_ratio

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.nozzle_mesh = Ring(height=outer_size[2], outer_top_radius=outer_size[0], inner_top_radius=inner_size[0],
                                outer_bottom_radius=outer_size[1], inner_bottom_radius=inner_size[1], x_z_ratio=x_z_ratio[0],
                                rotation=[0, np.pi / 2, 0])
        vertices_list.append(self.nozzle_mesh.vertices)
        faces_list.append(self.nozzle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.nozzle_mesh.vertices)

        nozzle_mesh_position = [0, inner_size[2] / 2, 0]

        self.nozzle_mesh = Cylinder(height=outer_size[2] - inner_size[2],
                                    top_radius=outer_size[0], bottom_radius=outer_size[0],
                                    top_radius_z=outer_size[0] * x_z_ratio[0], bottom_radius_z=outer_size[0] * x_z_ratio[0],
                                    position=nozzle_mesh_position)
        vertices_list.append(self.nozzle_mesh.vertices)
        faces_list.append(self.nozzle_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.nozzle_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Cap'


class Regular_cap(ConceptTemplate):
    def __init__(self, size, separation, position=[0, 0, 0], rotation=[0, 0, 0]):
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.separation = separation
        self.rotation = rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        self.mesh_1 = Cylinder(size[1], size[0])
        vertices_list.append(self.mesh_1.vertices)
        faces_list.append(self.mesh_1.faces + total_num_vertices)
        total_num_vertices += len(self.mesh_1.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order="YXZ", offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Cap'


def getvalue(self):
    attribute = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith("__")]
    for attr in attribute:
        print(f"{attr}: {getattr(self, attr)}")
