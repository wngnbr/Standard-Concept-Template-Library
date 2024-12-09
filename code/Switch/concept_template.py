import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh


class Standard_Base(ConceptTemplate):
    def __init__(self, size, has_back_part, back_part_size, back_part_offset, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.has_back_part = has_back_part
        self.size = size
        self.back_part_size = back_part_size
        self.back_part_offset = back_part_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        back_mesh_position = [0, 0, -size[2] / 2]
        self.back_mesh = Cuboid(size[1], size[0], size[2],
                                position=back_mesh_position)
        vertices_list.append(self.back_mesh.vertices)
        faces_list.append(self.back_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.back_mesh.vertices)

        if has_back_part[0]:
            back_mesh_position = [
                back_part_offset[0], 
                back_part_offset[1], 
                -size[2] - back_part_size[2] / 2
            ]
            self.back_mesh = Cuboid(back_part_size[1], back_part_size[0], back_part_size[2],
                                    position=back_mesh_position)
            vertices_list.append(self.back_mesh.vertices)
            faces_list.append(self.back_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.back_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Frame_Base(ConceptTemplate):
    def __init__(self, size, inner_size, inner_outer_offset, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.inner_size = inner_size
        self.inner_outer_offset = inner_outer_offset

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        back_mesh_rotation = [np.pi / 2, 0, 0]
        self.back_mesh = Rectangular_Ring(size[2], size[0], size[1],
                                          inner_size[0], inner_size[1],
                                          [inner_outer_offset[0], -inner_outer_offset[1]],
                                          rotation=back_mesh_rotation)
        vertices_list.append(self.back_mesh.vertices)
        faces_list.append(self.back_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.back_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Parallel_Base(ConceptTemplate):
    def __init__(self, size, sub_offset, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.sub_offset = sub_offset
        self.size = size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        left_mesh_position = [
            -size[0] / 2, 
            0, 
            -size[4] / 2
        ]
        self.left_mesh = Cuboid(size[2], size[0], size[4],
                                position=left_mesh_position)
        vertices_list.append(self.left_mesh.vertices)
        faces_list.append(self.left_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.left_mesh.vertices)

        right_mesh_position = [
            size[1] / 2, 
            sub_offset[0], 
            -size[4] / 2
        ]
        self.right_mesh = Cuboid(size[3], size[1], size[4],
                                 position=right_mesh_position)
        vertices_list.append(self.right_mesh.vertices)
        faces_list.append(self.right_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.right_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Body'


class Standard_Knob(ConceptTemplate):
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

        base_mesh_position = [0, 0, size[1] / 2]
        base_mesh_rotation = [np.pi / 2, 0, 0]
        self.base_mesh = Cylinder(size[1], size[0], size[0],
                                  position=base_mesh_position,
                                  rotation=base_mesh_rotation)
        vertices_list.append(self.base_mesh.vertices)
        faces_list.append(self.base_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.base_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, offset_first=True)
        
        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Knob'


class Round_Switch(ConceptTemplate):
    def __init__(self, number_of_switch, size, offset_1, offset_2, offset_3, offset_4, offset_Z, switch_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        switch_rotation = [x / 180 * np.pi for x in switch_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_switch = number_of_switch
        self.size = size
        self.offset_1 = offset_1
        self.offset_2 = offset_2
        self.offset_3 = offset_3
        self.offset_4 = offset_4
        self.offset_Z = offset_Z
        self.switch_rotation = switch_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(number_of_switch[0]):
            base_mesh_position = [
                locals()['offset_%d'%(i+1)][0],
                locals()['offset_%d'%(i+1)][1],
                offset_Z[0]
            ]
            base_mesh_rotation = [np.pi / 2 + switch_rotation[0], 0, 0]
            self.base_mesh = Cylinder(size[1], size[0], size[0],
                                      position=base_mesh_position,
                                      rotation=base_mesh_rotation)
            vertices_list.append(self.base_mesh.vertices)
            faces_list.append(self.base_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.base_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Switch'


class FlipX_Switch(ConceptTemplate):
    def __init__(self, number_of_switch, size, switch_rotation, separation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        switch_rotation = [x / 180 * np.pi for x in switch_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_switch = number_of_switch
        self.size = size
        self.separation = separation
        self.switch_rotation = switch_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(number_of_switch[0]):
            mesh_position = [
                (separation[0] + size[0]) * i, 
                0, 
                0
            ]
            mesh_rotation = [switch_rotation[0], 0, 0]
            self.mesh = Cuboid(size[1], size[0], size[2],
                               position=mesh_position,
                               rotation=mesh_rotation)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Switch'


class FlipY_Switch(ConceptTemplate):
    def __init__(self, number_of_switch, size, switch_rotation, separation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        switch_rotation = [x / 180 * np.pi for x in switch_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_switch = number_of_switch
        self.size = size
        self.separation = separation
        self.switch_rotation = switch_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(number_of_switch[0]):
            mesh_position = [
                (separation[0] + size[0]) * i, 
                0, 
                0
            ]
            mesh_rotation = [0, switch_rotation[0], 0]
            self.mesh = Cuboid(size[1], size[0], size[2],
                               position=mesh_position,
                               rotation=mesh_rotation)
            vertices_list.append(self.mesh.vertices)
            faces_list.append(self.mesh.faces + total_num_vertices)
            total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Switch'


class Lever_Switch(ConceptTemplate):
    def __init__(self, number_of_switch, base_size, main_size, inter_offset, switch_rotation, separation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        switch_rotation = [x / 180 * np.pi for x in switch_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_switch = number_of_switch
        self.base_size = base_size
        self.main_size = main_size
        self.inter_offset = inter_offset
        self.switch_rotation = switch_rotation
        self.separation = separation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for i in range(number_of_switch[0]):
            base_mesh_position = [
                (separation[0] + main_size[0]) * i,
                0,
                base_size[1] / 2,
            ]
            base_mesh_rotation = [np.pi / 2, 0, 0]
            self.base_mesh = Cylinder(base_size[1], base_size[0], base_size[0],
                                      position=base_mesh_position,
                                      rotation=base_mesh_rotation)
            vertices_list.append(self.base_mesh.vertices)
            faces_list.append(self.base_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.base_mesh.vertices)

            main_mesh_position = [
                (separation[0] + main_size[0]) * i,
                0,
                main_size[2] / 2,
            ]
            main_mesh_rotation = [np.pi / 2, 0, 0]
            self.main_mesh = Cylinder(main_size[2], main_size[1], main_size[0],
                                      position=main_mesh_position,
                                      rotation=main_mesh_rotation)
            self.main_mesh.vertices = apply_transformation(
                self.main_mesh.vertices,
                position=[0, 0, 0],
                rotation=[switch_rotation[0], 0, 0],
            )
            self.main_mesh.vertices = apply_transformation(
                self.main_mesh.vertices,
                position=[
                    inter_offset[0],
                    inter_offset[1],
                    inter_offset[2] + base_size[1],
                ],
                rotation=[0, 0, 0],
            )
            vertices_list.append(self.main_mesh.vertices)
            faces_list.append(self.main_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.main_mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Switch'


class Cuboidal_Plug(ConceptTemplate):
    def __init__(self, column_of_contact, row_of_contact, size, interval, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.column_of_contact = column_of_contact
        self.row_of_contact = row_of_contact
        self.size = size
        self.interval = interval

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for j in range(row_of_contact[0]):
            for i in range(column_of_contact[0]):
                mesh_position = [
                    (interval[0] + size[0]) * i,
                    (interval[1] + size[1]) * j,
                    -size[2] / 2,
                ],
                self.mesh = Cuboid(size[1], size[0], size[2],
                                   position=mesh_position)
                vertices_list.append(self.mesh.vertices)
                faces_list.append(self.mesh.faces + total_num_vertices)
                total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Plug'


class Standard_Plug(ConceptTemplate):
    def __init__(self, size, sub_offset, plug_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        plug_rotation = [x / 180 * np.pi for x in plug_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.sub_offset = sub_offset
        self.plug_rotation = plug_rotation
        self.size = size

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        mesh_position = [
            0, 
            0, 
            -size[2] / 2
        ]
        mesh_rotation = [0, 0, plug_rotation[0]]
        self.mesh = Cuboid(size[1], size[0], size[2],
                           position=mesh_position,
                           rotation=mesh_rotation)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        mesh_position = [
            - sub_offset[0] * np.cos(plug_rotation[0]) + sub_offset[1] * np.sin(-plug_rotation[0]),
            sub_offset[1] * np.cos(-plug_rotation[0]) - sub_offset[0] * np.sin(plug_rotation[0]),
            -size[2] / 2,
        ]
        mesh_rotation = [0, 0, plug_rotation[0] + plug_rotation[1]]
        self.mesh = Cuboid(size[1], size[0], size[2],
                           position=mesh_position,
                           rotation=mesh_rotation)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        mesh_position = [
            sub_offset[0] * np.cos(plug_rotation[0]) + sub_offset[1] * np.sin(-plug_rotation[0]),
            sub_offset[1] * np.cos(-plug_rotation[0]) + sub_offset[0] * np.sin(plug_rotation[0]),
            -size[2] / 2,
        ]
        mesh_rotation = [0, 0, plug_rotation[0] - plug_rotation[1]]
        self.mesh = Cuboid(size[1], size[0], size[2],
                           position=mesh_position,
                           rotation=mesh_rotation)
        vertices_list.append(self.mesh.vertices)
        faces_list.append(self.mesh.faces + total_num_vertices)
        total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Plug'


class Cylindrical_Plug(ConceptTemplate):
    def __init__(self, column_of_contact, row_of_contact, size, interval, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.column_of_contact = column_of_contact
        self.row_of_contact = row_of_contact
        self.size = size
        self.interval = interval

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0

        for j in range(row_of_contact[0]):
            for i in range(column_of_contact[0]):
                mesh_position = [
                    (interval[0] + size[0] * 2) * i,
                    (interval[1] + size[0] * 2) * j,
                    -size[1] / 2,
                ]
                mesh_rotation=[np.pi / 2, 0, 0]
                self.mesh = Cylinder(size[1], size[0], size[0],
                                     position=mesh_position,
                                     rotation=mesh_rotation)
                vertices_list.append(self.mesh.vertices)
                faces_list.append(self.mesh.faces + total_num_vertices)
                total_num_vertices += len(self.mesh.vertices)

        self.vertices = np.concatenate(vertices_list)
        self.faces = np.concatenate(faces_list)

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, offset_first=True)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Plug'