import numpy as np
from base_template import ConceptTemplate
from geometry_template import *
from utils import apply_transformation
from knowledge_utils import *
import trimesh


class Regular_Base(ConceptTemplate):
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
        
        self.back_mesh = Cuboid(size[1], size[0], size[2])
        vertices_list.append(self.back_mesh.vertices)
        faces_list.append(self.back_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.back_mesh.vertices)

        self.vertices=np.concatenate(vertices_list)
        self.faces=np.concatenate(faces_list)

        # Global Transformation
        self.vertices=apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Base'


class Regular_Screen(ConceptTemplate):
    def __init__(self, size, offset, screen_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):
        
        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        screen_rotation = [x / 180 * np.pi for x in screen_rotation]
        super().__init__(position, rotation)

        # Record Parameters
        self.size = size
        self.offset = offset
        self.screen_rotation = screen_rotation

        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0
        
        back_mesh_position = [
            0, 
            offset[0] + size[1] * np.cos(screen_rotation[0]) / 2,
            offset[1]
        ]
        back_mesh_rotation = [screen_rotation[0], 0, 0]
        self.back_mesh = Cuboid(size[1], size[0], size[2],
                                position=back_mesh_position,
                                rotation=back_mesh_rotation)
        vertices_list.append(self.back_mesh.vertices)
        faces_list.append(self.back_mesh.faces + total_num_vertices)
        total_num_vertices += len(self.back_mesh.vertices)
        
        self.vertices=np.concatenate(vertices_list)
        self.faces=np.concatenate(faces_list)

        # Global Transformation
        self.vertices=apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Screen'
        

class Cuboidal_Connector(ConceptTemplate):
    def __init__(self, number_of_connector, size, separation,offset, connector_rotation, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        connector_rotation = [x / 180 * np.pi for x in connector_rotation]    
        super().__init__(position, rotation)

        # Record Parameters
        self.number_of_connector = number_of_connector
        self.size = size
        self.offset = offset
        self.separation = separation
        self.connector_rotation = connector_rotation
        
        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0
        
        for i in range(number_of_connector[0]):
            back_mesh_position = [
                offset[0] + i * size[0] + sum(separation[0:i]) + size[0] / 2, 
                offset[1], 
                offset[2]
            ]
            back_mesh_rotation = [connector_rotation[0], 0, 0]
            self.back_mesh = Cuboid(size[1], size[0], size[2],
                                    position=back_mesh_position,
                                    rotation=back_mesh_rotation)
            vertices_list.append(self.back_mesh.vertices)
            faces_list.append(self.back_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.back_mesh.vertices)

        self.vertices=np.concatenate(vertices_list)
        self.faces=np.concatenate(faces_list)

        # Global Transformation
        self.vertices=apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Connector'

class Cylindrical_Connector(ConceptTemplate):
    def __init__(self, number_of_connector, size, separation, offset, position=[0, 0, 0], rotation=[0, 0, 0]):

        # Process rotation param
        rotation = [x / 180 * np.pi for x in rotation]
        super().__init__(position, rotation)  

        # Record Parameters
        self.number_of_connector = number_of_connector
        self.size = size
        self.offset = offset
        self.separation = separation
        
        # Instantiate component geometries
        vertices_list = []
        faces_list = []
        total_num_vertices = 0
         
        for i in range(number_of_connector[0]):
            back_mesh_position = [
                offset[0] + i * size[1] + sum(separation[0:i]) + size[1] / 2, 
                offset[1], 
                offset[2]
            ]
            back_mesh_rotation = [0, 0, np.pi / 2]
            self.back_mesh = Cylinder(size[1], size[0], size[0],
                                      position=back_mesh_position,
                                      rotation=back_mesh_rotation)
            vertices_list.append(self.back_mesh.vertices)
            faces_list.append(self.back_mesh.faces + total_num_vertices)
            total_num_vertices += len(self.back_mesh.vertices)

        self.vertices=np.concatenate(vertices_list)
        self.faces=np.concatenate(faces_list)
        
        # Global Transformation
        self.vertices=apply_transformation(self.vertices, position, rotation)

        self.overall_obj_mesh = trimesh.Trimesh(self.vertices, self.faces)
        self.overall_obj_pts = np.array(self.overall_obj_mesh.sample(SAMPLENUM))

        self.semantic = 'Connector'
        

