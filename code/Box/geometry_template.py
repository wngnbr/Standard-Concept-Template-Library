import numpy as np
from base_template import GeometryTemplate
from utils import apply_transformation

class Cuboid(GeometryTemplate):
    def __init__(self, height, top_length, top_width = None, bottom_length = None, bottom_width = None, top_offset = [0, 0], back_height = None, position = [0, 0, 0], rotation = [0, 0, 0], rotation_order = "XYZ"):
        """
        :param height: height of the front surface of the cuboid in the Y-axis direction
        :param top_length: length of the top surface of the cuboid in the X-axis direction
        :param top_width: width of the top surface of the cuboid in the Z-axis direction
        :param bottom_length: length of the bottom surface of the cuboid in the X-axis direction
        :param bottom_width: width of the bottom surface of the cuboid in the Z-axis direction
        :param top_offset: offset between the upper and lower surface of the cuboid in the X-axis and Z-axis directions
        :param back_height: height of the back surface of the cuboid in the Y-axis direction
        :param position: position (x, y, z) of the cuboid
        :param rotation: rotation of the cuboid, represented via Euler angles (x, y, z)
        :param rotation_order: rotation order of the three rotation axes of the cuboid

        """

        # Filling Missing Values
        if top_width == None:
            top_width = top_length
        if bottom_length == None:
            bottom_length = top_length
        if bottom_width == None:
            bottom_width = top_width
        if back_height == None:
            back_height = height

        super().__init__(position, rotation, rotation_order)

        # Record Parameters
        self.height = height
        self.top_length = top_length
        self.top_width = top_width
        self.bottom_length = bottom_length
        self.bottom_width = bottom_width
        self.top_offset = top_offset
        self.back_height = back_height
        self.position = position
        self.rotation = rotation
        self.rotation_order = rotation_order
            
        # Manually Defined Default Template Instance 
        self.vertices = np.array([
            [-1 / 2, 1 / 2, 1 / 2],
            [1 / 2, 1 / 2, 1 / 2],
            [-1 / 2, 1 / 2, -1 / 2],
            [1 / 2, 1 / 2, -1 / 2],
            [-1 / 2, -1 / 2, 1 / 2],
            [1 / 2, -1 / 2, 1 / 2],
            [-1 / 2, -1 / 2, -1 / 2],
            [1 / 2, -1 / 2, -1 / 2]
        ])

        self.faces = np.array([
            [0, 1, 2], [1, 3, 2],
            [4, 6, 5], [5, 6, 7],
            [0, 4, 5], [0, 5, 1],
            [2, 7, 6], [2, 3, 7],
            [0, 6, 4], [0, 2, 6],
            [1, 5, 7], [1, 7, 3]
        ])

        # Differentiable Deformation
        vertices_resize = np.array([
            [top_length, height, top_width], 
            [top_length, height, top_width], 
            [top_length, back_height, top_width],
            [top_length, back_height, top_width],
            [bottom_length, height, bottom_width],
            [bottom_length, height, bottom_width],
            [bottom_length, back_height, bottom_width],
            [bottom_length, back_height, bottom_width]
        ])
        self.vertices = self.vertices * vertices_resize

        vertices_offset = np.array([
            [top_offset[0], 0, top_offset[1]], 
            [top_offset[0], 0, top_offset[1]], 
            [top_offset[0], (back_height - height) / 2, top_offset[1]],
            [top_offset[0], (back_height - height) / 2, top_offset[1]],
            [0, 0, 0],
            [0, 0, 0],
            [0, (back_height - height) / 2, 0],
            [0, (back_height - height) / 2, 0]
        ])
        self.vertices = self.vertices + vertices_offset

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order)


class Sphere(GeometryTemplate):
    def __init__(self, radius, top_angle = 0, bottom_angle = np.pi, radius_y = None, radius_z = None, longitude_angle = np.pi*2, position = [0, 0, 0], rotation = [0, 0, 0], rotation_order = "XYZ"):
        """
        :param radius: length of the half-axis of the sphere in the X-axis direction
        :param top_angle: latitude starting angle of the sphere
        :param bottom_angle: latitude ending angle of the sphere
        :param radius_y: length of the half-axis of the sphere in the Y-axis direction
        :param radius_z: length of the half-axis of the sphere in the Z-axis direction
        :param longitude_angle: longitude covered angle of the sphere
        :param position: position (x, y, z) of the sphere
        :param rotation: rotation of the sphere, represented via Euler angles (x, y, z)
        :param rotation_order: rotation order of the three rotation axes of the sphere

        """

        # Filling Missing Values
        if radius_y == None:
            radius_y = radius
        if radius_z == None:
            radius_z = radius

        super().__init__(position, rotation, rotation_order)

        # Record Parameters
        self.radius = radius
        self.top_angle = top_angle
        self.bottom_angle = bottom_angle
        self.radius_y = radius_y
        self.radius_z = radius_z
        self.longitude_angle = longitude_angle
        self.position = position
        self.rotation = rotation
        self.rotation_order = rotation_order
            
        # Manually Defined Default Template Instance 
        vertices = []
        num_of_segment_latitude = 64
        num_of_segment_longitude = 64
        for i in range(num_of_segment_latitude+1):
            for j in range(num_of_segment_longitude+1):
                vertices.append([0, 0, 0])
        vertices.append([0, 0, 0])
        vertices.append([0, 0, 0])
        self.vertices = np.array(vertices)

        faces = []
        for i in range(num_of_segment_latitude):
            for j in range(num_of_segment_longitude):
                faces.append([(num_of_segment_longitude+1)* i + j, (num_of_segment_longitude+1) * i + j + 1, (num_of_segment_longitude+1) * (i + 1) + j + 1])
                faces.append([(num_of_segment_longitude+1)* i + j, (num_of_segment_longitude+1) * (i + 1) + j + 1, (num_of_segment_longitude+1) * (i + 1) + j])

        for i in range(num_of_segment_longitude):
            faces.append([(num_of_segment_latitude+1)*(num_of_segment_longitude+1), i + 1, i])
            faces.append([(num_of_segment_latitude+1)*(num_of_segment_longitude+1) + 1, num_of_segment_latitude * (num_of_segment_longitude+1) + i, num_of_segment_latitude * (num_of_segment_longitude+1) + i + 1])
        self.faces = np.array(faces)

        # Differentiable Deformation
        radius_offset = None
        for i in range(num_of_segment_latitude+1):
            rotation_radius_tmp = (bottom_angle - top_angle) / num_of_segment_latitude * i + top_angle
            radius_offset_tmp = np.array([radius * np.sin(rotation_radius_tmp), radius_y * np.cos(rotation_radius_tmp), radius_z * np.sin(rotation_radius_tmp)])
            radius_offset_tmp = np.tile(radius_offset_tmp[None, :], (num_of_segment_longitude+1, 1))
            if radius_offset is None:
                radius_offset = radius_offset_tmp
            else:
                radius_offset = np.concatenate((radius_offset, radius_offset_tmp), axis=0)

        central_angle = None
        for i in range(num_of_segment_longitude+1):
            rotation_longitude_tmp = longitude_angle / num_of_segment_longitude * i
            central_angle_tmp = np.array([np.cos(rotation_longitude_tmp), 1, np.sin(rotation_longitude_tmp)])
            if central_angle is None:
                central_angle = central_angle_tmp
            elif i == 1:
                central_angle = np.array([central_angle, central_angle_tmp])
            else:
                central_angle = np.concatenate((central_angle, central_angle_tmp[None, :]), axis=0)
        central_angle = np.tile(central_angle, (int((self.vertices.shape[0]) / ((num_of_segment_longitude+1))), 1))
        radius_offset *= central_angle

        top_center = np.array([0, radius_y * np.cos(top_angle), 0])
        bottom_center = np.array([0, radius_y * np.cos(bottom_angle), 0])
        center = np.array([top_center, bottom_center])
        radius_offset = np.concatenate((radius_offset, center), axis=0)

        self.vertices = self.vertices + radius_offset

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order)

class Cylinder(GeometryTemplate):
    def __init__(self, height, top_radius, bottom_radius = None, top_radius_z = None, bottom_radius_z = None, is_half = False, is_quarter = False, position = [0, 0, 0], rotation = [0, 0, 0], rotation_order = "XYZ"):
        """
        :param height: height of the cylinder in the Y-axis direction
        :param top_radius: radius of the top surface of the cylinder in the X-axis direction
        :param bottom_radius: radius of the bottom surface of the cylinder in the X-axis direction
        :param top_radius_z: radius of the top surface of the cylinder in the Z-axis direction
        :param bottom_radius_z: radius of the bottom surface of the cylinder in the Z-axis direction
        :param is_half: whether the cylinder is half
        :param is_quarter: whether the cylinder is quarter
        :param position: position (x, y, z) of the cylinder
        :param rotation: rotation of the cylinder, represented via Euler angles (x, y, z)
        :param rotation_order: rotation order of the three rotation axes of the cylinder

        """

        # Filling Missing Values
        if bottom_radius == None:
            bottom_radius = top_radius
        if top_radius_z == None:
            top_radius_z = top_radius
        if bottom_radius_z == None:
            bottom_radius_z = bottom_radius

        super().__init__(position, rotation, rotation_order)

        # Record Parameters
        self.height = height
        self.top_radius = top_radius
        self.bottom_radius = bottom_radius
        self.top_radius_z = top_radius_z
        self.bottom_radius_z = bottom_radius_z
        self.is_half = is_half
        self.is_quarter = is_quarter
        self.position = position
        self.rotation = rotation
        self.rotation_order = rotation_order
            
        # Manually Defined Default Template Instance 
        vertices = []
        vertices.append([0, 1 / 2, 0])
        vertices.append([0, -1 / 2, 0])
        num_of_segment = 256
        for i in range(num_of_segment+1):
            rotation_tmp = np.pi * 2 / num_of_segment * i
            if is_half:
                rotation_tmp = rotation_tmp / 2
            elif is_quarter:
                rotation_tmp = rotation_tmp / 4
            vertices.append([np.cos(rotation_tmp), 1 / 2, np.sin(rotation_tmp)])
            vertices.append([np.cos(rotation_tmp), -1 / 2, np.sin(rotation_tmp)])
        self.vertices = np.array(vertices)

        faces = []
        faces.append([1, 0, 3])
        faces.append([0, 2, 3])
        for i in range(num_of_segment):
            faces.append([2*i+2, 2*i+4, 2*i+3])
            faces.append([2*i+5, 2*i+3, 2*i+4])
            faces.append([2*i+2, 0, 2*i+4])
            faces.append([1, 2*i+3, 2*i+5])
        faces.append([0, 1, 2*(num_of_segment+1)+1])
        faces.append([2*(num_of_segment+1)+1, 0, 2*(num_of_segment+1)])
        self.faces = np.array(faces)

        # Differentiable Deformation
        vertices_resize = np.array([
            [top_radius, height, top_radius_z], 
            [bottom_radius, height, bottom_radius_z]
        ])
        vertices_resize = np.tile(vertices_resize, (int(self.vertices.shape[0] / 2), 1))
        self.vertices = self.vertices * vertices_resize

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order)


class Trianguler_Prism(GeometryTemplate):
    def __init__(self, height, top_radius, bottom_radius = None, position = [0, 0, 0], rotation = [0, 0, 0], rotation_order = "XYZ"):
        """
        :param height: height of the trianguler in the Y-axis direction
        :param top_radius: length from the base vertex to the center of the top surface of the trianguler_prism
        :param bottom_radius: length from the base vertex to the center of the bottom surface of the trianguler_prism
        :param position: position (x, y, z) of the trianguler
        :param rotation: rotation of the trianguler, represented via Euler angles (x, y, z)
        :param rotation_order: rotation order of the three rotation axes of the trianguler

        """

        # Filling Missing Values
        if bottom_radius == None:
            bottom_radius = top_radius

        super().__init__(position, rotation, rotation_order)

        # Record Parameters
        self.height = height
        self.top_radius = top_radius
        self.bottom_radius = bottom_radius
        self.position = position
        self.rotation = rotation
        self.rotation_order = rotation_order
            
        # Manually Defined Default Template Instance 
        vertices = []
        vertices.append([0, 1 / 2, 0])
        vertices.append([0, -1 / 2, 0])
        num_of_segment = 3
        for i in range(num_of_segment+1):
            rotation_tmp = np.pi * 2 / num_of_segment * i
            vertices.append([np.cos(rotation_tmp), 1 / 2, np.sin(rotation_tmp)])
            vertices.append([np.cos(rotation_tmp), -1 / 2, np.sin(rotation_tmp)])
        self.vertices = np.array(vertices)

        faces = []
        for i in range(num_of_segment):
            faces.append([2*i+2, 2*i+4, 2*i+3])
            faces.append([2*i+5, 2*i+3, 2*i+4])
            faces.append([2*i+2, 0, 2*i+4])
            faces.append([1, 2*i+3, 2*i+5])
        self.faces = np.array(faces)

        # Differentiable Deformation
        vertices_resize = np.array([
            [top_radius, height, top_radius], 
            [bottom_radius, height, bottom_radius]
        ])
        vertices_resize = np.tile(vertices_resize, (int(self.vertices.shape[0] / 2), 1))
        self.vertices = self.vertices * vertices_resize

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order)


class Cone(GeometryTemplate):
    def __init__(self, radius, height, tip_offset = [0, 0], radius_z = None, position = [0, 0, 0], rotation = [0, 0, 0], rotation_order = "XYZ"):
        """
        :param radius: radius of the bottom surface of the cone in the X-axis direction
        :param height: height of the cone in the Y-axis direction
        :param tip_offset: offset between the apex of the cone and the center of the bottom circle in the X-axis and Z-axis direction
        :param radius_z: radius of the bottom surface of the cone in the Z-axis direction
        :param position: position (x, y, z) of the cone
        :param rotation: rotation of the cone, represented via Euler angles (x, y, z)
        :param rotation_order: rotation order of the three rotation axes of the cone

        """

        # Filling Missing Values
        if radius_z == None:
            radius_z = radius

        super().__init__(position, rotation, rotation_order)

        # Record Parameters
        self.radius = radius
        self.height = height
        self.tip_offset = tip_offset
        self.radius_z = radius_z
        self.position = position
        self.rotation = rotation
        self.rotation_order = rotation_order
            
        # Manually Defined Default Template Instance 
        vertices = []
        num_of_segment = 256
        vertices.append([0, 0, 0])
        vertices.append([0, 0, 0])
        for i in range(num_of_segment+1):
            vertices.append([0, 0, 0])
        self.vertices = np.array(vertices)

        faces = []
        for i in range(num_of_segment):
            faces.append([0, i + 2, i + 3])
            faces.append([1, i + 3, i + 2])
        self.faces = np.array(faces)

        # Differentiable Deformation
        vertices_resize = np.array([
            [0, 0, 0], 
            [tip_offset[0], height, tip_offset[1]]
        ])
        for i in range(num_of_segment+1):
            rotation_tmp = np.pi * 2 / num_of_segment * i
            size_tmp = np.array([radius * np.cos(rotation_tmp), 0, radius_z * np.sin(rotation_tmp)])
            vertices_resize = np.concatenate([vertices_resize, size_tmp[None, :]], axis=0)
        self.vertices = self.vertices + vertices_resize

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order)


class Rectangular_Ring(GeometryTemplate):
    def __init__(self, front_height, outer_top_length, outer_top_width, inner_top_length, inner_top_width, inner_offset = [0, 0], outer_bottom_length = None, outer_bottom_width = None, inner_bottom_length = None, inner_bottom_width = None, back_height = None, top_bottom_offset = [0, 0], position = [0, 0, 0], rotation = [0, 0, 0], rotation_order = "XYZ"):
        """
        :param front_height: height of the outer part of the rectangular ring in the Y-axis direction
        :param outer_top_length: length of the top surface of the outer part of the rectangular ring in the X-axis direction
        :param outer_top_width: width of the top surface of the outer part of the rectangular ring in the Z-axis direction
        :param inner_top_length: length of the top surface of the inner part of the rectangular ring in the X-axis direction
        :param inner_top_width: width of the top surface of the inner part of the rectangular ring in the Z-axis direction
        :param inner_offset: offset between the outer and inner part of the rectangular ring in the X-axis and Z-axis directions
        :param outer_bottom_length: length of the bottom surface of the outer part of the rectangular ring in the X-axis direction
        :param outer_bottom_width: width of the bottom surface of the outer part of the rectangular ring in the Z-axis direction
        :param inner_bottom_length: length of the bottom surface of the inner part of the rectangular ring in the X-axis direction
        :param inner_bottom_width: width of the bottom surface of the inner part of the rectangular ring in the Z-axis direction
        :param back_height: height of the back surface of the outer part of the rectangular ring in the Y-axis direction
        :param top_bottom_offset: offset between the upper and lower surface of the outer part of the rectangular ring in the X-axis and Z-axis directions
        :param position: position (x, y, z) of the rectangular ring
        :param rotation: rotation of the rectangular ring, represented via Euler angles (x, y, z)
        :param rotation_order: rotation order of the three rotation axes of the rectangular ring

        """

        # Filling Missing Values
        if outer_bottom_length == None:
            outer_bottom_length = outer_top_length
        if outer_bottom_width == None:
            outer_bottom_width = outer_top_width
        if inner_bottom_length == None:
            inner_bottom_length = inner_top_length
        if inner_bottom_width == None:
            inner_bottom_width = inner_top_width
        if back_height == None:
            back_height = front_height

        super().__init__(position, rotation, rotation_order)

        # Record Parameters
        self.front_height = front_height
        self.outer_top_length = outer_top_length
        self.outer_top_width = outer_top_width
        self.inner_top_length = inner_top_length
        self.inner_top_width = inner_top_width
        self.inner_offset = inner_offset
        self.outer_bottom_length = outer_bottom_length
        self.outer_bottom_width = outer_bottom_width
        self.inner_bottom_length = inner_bottom_length
        self.inner_bottom_width = inner_bottom_width
        self.back_height = back_height
        self.top_bottom_offset = top_bottom_offset
        self.position = position
        self.rotation = rotation
        self.rotation_order = rotation_order
            
        # Manually Defined Default Template Instance 
        self.vertices = np.array([
            [-1 / 2, 1 / 2, 1 / 2], # outer_box
            [1 / 2, 1 / 2, 1 / 2],
            [-1 / 2, 1 / 2, -1 / 2],
            [1 / 2, 1 / 2, -1 / 2],
            [-1 / 2, -1 / 2, 1 / 2],
            [1 / 2, -1 / 2, 1 / 2],
            [-1 / 2, -1 / 2, -1 / 2],
            [1 / 2, -1 / 2, -1 / 2],
            [-1 / 2, 1 / 2, 1 / 2], # inner_box
            [1 / 2, 1 / 2, 1 / 2],
            [-1 / 2, 1 / 2, -1 / 2],
            [1 / 2, 1 / 2, -1 / 2],
            [-1 / 2, -1 / 2, 1 / 2],
            [1 / 2, -1 / 2, 1 / 2],
            [-1 / 2, -1 / 2, -1 / 2],
            [1 / 2, -1 / 2, -1 / 2]
        ])

        self.faces = np.array([
            [0, 4, 5], [0, 5, 1], # outer_box
            [2, 7, 6], [2, 3, 7],
            [0, 6, 4], [0, 2, 6],
            [1, 5, 7], [1, 7, 3],
            [4+8, 0+8, 5+8], [5+8, 0+8, 1+8], # inner_box
            [7+8, 2+8, 6+8], [3+8, 2+8, 7+8],
            [6+8, 0+8, 4+8], [2+8, 0+8, 6+8],
            [5+8, 1+8, 7+8], [7+8, 1+8, 3+8],
            [0, 1, 1+8], [1+8, 0+8, 0], # top-outer-inner
            [1, 3, 3+8], [3+8, 1+8, 1],
            [3, 2, 2+8], [2+8, 3+8, 3],
            [2, 0, 0+8], [0+8, 2+8, 2],
            [5, 4, 5+8], [4+8, 5+8, 4], # bottom-outer-inner
            [7, 5, 7+8], [5+8, 7+8, 5],
            [6, 7, 6+8], [7+8, 6+8, 7],
            [4, 6, 4+8], [6+8, 4+8, 6],
        ])

        # Differentiable Deformation
        vertices_resize = np.array([
            [outer_top_length, front_height, outer_top_width], 
            [outer_top_length, front_height, outer_top_width], 
            [outer_top_length, back_height, outer_top_width],
            [outer_top_length, back_height, outer_top_width],
            [outer_bottom_length, front_height, outer_bottom_width],
            [outer_bottom_length, front_height, outer_bottom_width],
            [outer_bottom_length, back_height, outer_bottom_width],
            [outer_bottom_length, back_height, outer_bottom_width],
            [inner_top_length, front_height, inner_top_width], 
            [inner_top_length, front_height, inner_top_width], 
            [inner_top_length, back_height, inner_top_width],
            [inner_top_length, back_height, inner_top_width],
            [inner_bottom_length, front_height, inner_bottom_width],
            [inner_bottom_length, front_height, inner_bottom_width],
            [inner_bottom_length, back_height, inner_bottom_width],
            [inner_bottom_length, back_height, inner_bottom_width]
        ])
        self.vertices = self.vertices * vertices_resize

        vertices_offset = np.array([
            [0, 0, 0], 
            [0, 0, 0], 
            [0, (back_height - front_height) / 2, 0],
            [0, (back_height - front_height) / 2, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, (back_height - front_height) / 2, 0],
            [0, (back_height - front_height) / 2, 0],
            [inner_offset[0], 0, inner_offset[1]], 
            [inner_offset[0], 0, inner_offset[1]], 
            [inner_offset[0], (back_height - front_height) / 2, inner_offset[1]],
            [inner_offset[0], (back_height - front_height) / 2, inner_offset[1]],
            [inner_offset[0], 0, inner_offset[1]],
            [inner_offset[0], 0, inner_offset[1]],
            [inner_offset[0], (back_height - front_height) / 2, inner_offset[1]],
            [inner_offset[0], (back_height - front_height) / 2, inner_offset[1]]
        ])
        self.vertices = self.vertices + vertices_offset

        vertices_top_offset = np.array([
            [top_bottom_offset[0], 0, top_bottom_offset[1]], 
            [top_bottom_offset[0], 0, top_bottom_offset[1]],
            [top_bottom_offset[0], 0, top_bottom_offset[1]],
            [top_bottom_offset[0], 0, top_bottom_offset[1]],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [top_bottom_offset[0], 0, top_bottom_offset[1]], 
            [top_bottom_offset[0], 0, top_bottom_offset[1]],
            [top_bottom_offset[0], 0, top_bottom_offset[1]],
            [top_bottom_offset[0], 0, top_bottom_offset[1]],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ])
        self.vertices = self.vertices + vertices_top_offset

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order)


class Ring(GeometryTemplate):
    def __init__(self, height, outer_top_radius, inner_top_radius, exist_angle = np.pi*2, outer_bottom_radius = None, inner_bottom_radius = None, back_height = None, generatrix_offset = 0, x_z_ratio = 1, inner_x_z_ratio = None, inner_offset = [0, 0], position = [0, 0, 0], rotation = [0, 0, 0], rotation_order = "XYZ"):
        """
        :param height: height of the ring in the Y-axis direction
        :param outer_top_radius: radius of the top surface of the outer part of the ring in the X-axis direction
        :param inner_top_radius: radius of the top surface of the inner part of the ring in the X-axis direction
        :param exist_angle: covered angle of the ring
        :param outer_bottom_radius: radius of the bottom surface of the outer part of the ring in the X-axis direction
        :param inner_bottom_radius: radius of the bottom surface of the inner part of the ring in the X-axis direction
        :param back_height: height of the back generatrix of the outer part of the ring in the Y-axis direction
        :param generatrix_offset: offset between the front and back generatrix of the outer part of the ring in the Y-axis direction
        :param x_z_ratio: ratio of outer radius of the ring in the X-axis and Z-axis directions
        :param inner_x_z_ratio: ratio of inner radius of the ring in the X-axis and Z-axis directions
        :param inner_offset: offset between the outer and inner part of the ring in the Y-axis and Z-axis directions
        :param position: position (x, y, z) of the ring
        :param rotation: rotation of the ring, represented via Euler angles (x, y, z)
        :param rotation_order: rotation order of the three rotation axes of the ring

        """

        # Filling Missing Values
        if outer_bottom_radius == None:
            outer_bottom_radius = outer_top_radius
        if inner_bottom_radius == None:
            inner_bottom_radius = inner_top_radius
        if back_height == None:
            back_height = height
        if inner_x_z_ratio == None:
            inner_x_z_ratio = x_z_ratio

        super().__init__(position, rotation, rotation_order)

        # Record Parameters
        self.height = height
        self.outer_top_radius = outer_top_radius
        self.inner_top_radius = inner_top_radius
        self.exist_angle = exist_angle
        self.outer_bottom_radius = outer_bottom_radius
        self.inner_bottom_radius = inner_bottom_radius
        self.back_height = back_height
        self.generatrix_offset = generatrix_offset
        self.x_z_ratio = x_z_ratio
        self.inner_x_z_ratio = inner_x_z_ratio
        self.inner_offset = inner_offset
        self.position = position
        self.rotation = rotation
        self.rotation_order = rotation_order
            
        # Manually Defined Default Template Instance 
        vertices = []
        num_of_segment = 256
        for i in range(num_of_segment+1):
            vertices.append([0, 0, 0])
            vertices.append([0, 0, 0])
            vertices.append([0, 0, 0])
            vertices.append([0, 0, 0])
        self.vertices = np.array(vertices)

        faces = []
        faces.append([0, 1, 2])
        faces.append([1, 3, 2])
        for i in range(num_of_segment):
            faces.append([4 * i + 0, 4 * i + 2, 4 * (i + 1) + 2])
            faces.append([4 * i + 0, 4 * (i + 1) + 2, 4 * (i + 1)])
            faces.append([4 * i + 1, 4 * (i + 1) + 3, 4 * i + 3])
            faces.append([4 * i + 1, 4 * (i + 1) + 1, 4 * (i + 1) + 3])
            faces.append([4 * i + 1, 4 * i + 0, 4 * (i + 1) + 0])
            faces.append([4 * i + 1, 4 * (i + 1) + 0, 4 * (i + 1) + 1])
            faces.append([4 * (i + 1) + 2, 4 * i + 2, 4 * i + 3])
            faces.append([4 * (i + 1) + 2, 4 * i + 3, 4 * (i + 1) + 3])
        faces.append([num_of_segment*4 + 0, num_of_segment*4 + 2, num_of_segment*4 + 1])
        faces.append([num_of_segment*4 + 1, num_of_segment*4 + 2, num_of_segment*4 + 3])
        self.faces = np.array(faces)

        # Differentiable Deformation
        size = None
        for i in range(num_of_segment+1):
            rotation_tmp = exist_angle / num_of_segment * i
            top_y = height / 2 + (back_height / 2 + generatrix_offset - height / 2) / (outer_top_radius * 2) * (outer_top_radius * (1 - np.cos(rotation_tmp)))
            bottom_y = -height / 2 + (-back_height / 2 + generatrix_offset + height / 2) / (outer_bottom_radius * 2) * (outer_bottom_radius * (1 - np.cos(rotation_tmp)))
            top_outer_vert = np.array([outer_top_radius * x_z_ratio * np.cos(rotation_tmp), top_y, outer_top_radius * np.sin(rotation_tmp)])
            bottom_outer_vert = np.array([outer_bottom_radius * x_z_ratio * np.cos(rotation_tmp), bottom_y, outer_bottom_radius * np.sin(rotation_tmp)])

            sur_delta_y = top_y - height / 2
            sur_delta_y_new = min(sur_delta_y, (back_height / 2 + generatrix_offset - height / 2) - sur_delta_y)
            sur_y_gap = (back_height / 2 + generatrix_offset - height / 2) - sur_delta_y_new * 2
            section_chord = np.sqrt((outer_top_radius * 2) * (outer_top_radius * 2) + sur_y_gap * sur_y_gap)
            delta_chord = section_chord * (outer_top_radius - inner_top_radius) / (outer_top_radius * 2)
            delta_y = np.sqrt(np.abs(delta_chord * delta_chord - (outer_top_radius - inner_top_radius) * (outer_top_radius - inner_top_radius)))
            if sur_delta_y > ((back_height / 2 + generatrix_offset - height / 2) - sur_delta_y):
                delta_y = -delta_y
            top_inner_vert = np.array([inner_top_radius * inner_x_z_ratio * np.cos(rotation_tmp), top_y + delta_y, inner_top_radius * np.sin(rotation_tmp)])
            top_inner_vert += np.array([0, inner_offset[1], -inner_offset[0]])
            
            sur_delta_y = bottom_y + height / 2
            sur_delta_y_new = max(sur_delta_y, (-back_height / 2 + generatrix_offset + height / 2) - sur_delta_y)
            sur_y_gap = (-back_height / 2 + generatrix_offset + height / 2) - sur_delta_y_new * 2
            section_chord = np.sqrt((outer_bottom_radius * 2) * (outer_bottom_radius * 2) + sur_y_gap * sur_y_gap)
            delta_chord = section_chord * (outer_bottom_radius - inner_bottom_radius) / (outer_bottom_radius * 2)
            delta_y = np.sqrt(np.abs(delta_chord * delta_chord - (outer_bottom_radius - inner_bottom_radius) * (outer_bottom_radius - inner_bottom_radius)))
            if sur_delta_y > ((-back_height / 2 + generatrix_offset + height / 2) - sur_delta_y):
                delta_y = -delta_y
            bottom_inner_vert = np.array([inner_bottom_radius * inner_x_z_ratio * np.cos(rotation_tmp), bottom_y + delta_y, inner_bottom_radius * np.sin(rotation_tmp)])
            bottom_inner_vert += np.array([0, inner_offset[1], -inner_offset[0]])
            
            size_tmp = np.array([top_outer_vert, bottom_outer_vert, top_inner_vert, bottom_inner_vert])
            if size is None:
                size = size_tmp
            else:
                size = np.concatenate((size, size_tmp), axis=0)

        self.vertices = self.vertices + size

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order)


class Torus(GeometryTemplate):
    def __init__(self, central_radius, start_torus_radius, exist_angle = np.pi*2, end_torus_radius = None, position = [0, 0, 0], rotation = [0, 0, 0], rotation_order = "XYZ"):
        """
        :param central_radius: central radius of the torus
        :param start_torus_radius: starting torus radius of the torus
        :param exist_angle: covered angle of the torus
        :param end_torus_radius: ending torus radius of the torus
        :param position: position (x, y, z) of the torus
        :param rotation: rotation of the torus, represented via Euler angles (x, y, z)
        :param rotation_order: rotation order of the three rotation axes of the torus

        """

        # Filling Missing Values
        if end_torus_radius == None:
            end_torus_radius = start_torus_radius

        super().__init__(position, rotation, rotation_order)

        # Record Parameters
        self.central_radius = central_radius
        self.start_torus_radius = start_torus_radius
        self.exist_angle = exist_angle
        self.end_torus_radius = end_torus_radius
        self.position = position
        self.rotation = rotation
        self.rotation_order = rotation_order
            
        # Manually Defined Default Template Instance 
        vertices = []
        num_of_segment_center = 64
        num_of_segment_trous = 64
        vertices.append([0, 0, 0])
        vertices.append([0, 0, 0])
        for i in range(num_of_segment_center+1):
            for j in range(num_of_segment_trous+1):
                vertices.append([0, 0, 0])
        self.vertices = np.array(vertices)

        faces = []
        for i in range(num_of_segment_trous):
            faces.append([0, 2 + i + 1, 2 + i])
            faces.append([1, 2 + (num_of_segment_center+1) * num_of_segment_trous + i, 2 + (num_of_segment_center+1) * num_of_segment_trous + i + 1])

        for i in range(num_of_segment_center):
            for j in range(num_of_segment_trous):
                faces.append([2 + (num_of_segment_trous + 1) * i + j + 1, 2 + (num_of_segment_trous + 1) * (i + 1) + j, 2 + (num_of_segment_trous + 1) * i + j])
                faces.append([2 + (num_of_segment_trous + 1) * i + j + 1, 2 + (num_of_segment_trous + 1) * (i + 1) + j + 1, 2 + (num_of_segment_trous + 1) * (i + 1) + j])
        self.faces = np.array(faces)

        # Differentiable Deformation
        size = np.array([
            [central_radius, 0, 0],
            [central_radius * np.cos(exist_angle), 0, central_radius * np.sin(exist_angle)]
        ])
        for i in range(num_of_segment_center+1):
            rotation_center_tmp = exist_angle / num_of_segment_center * i
            outer_trous_radius_tmp = (start_torus_radius * (num_of_segment_center - i) + end_torus_radius * i) / num_of_segment_center
            for j in range(num_of_segment_trous+1):
                rotation_trous_tmp = np.pi * 2 / num_of_segment_trous * j
                outer_total_length = central_radius + outer_trous_radius_tmp * np.cos(rotation_trous_tmp)
                outer_size_tmp = np.array([outer_total_length * np.cos(rotation_center_tmp), outer_trous_radius_tmp * np.sin(rotation_trous_tmp), outer_total_length * np.sin(rotation_center_tmp)])
                size_tmp = outer_size_tmp[None, :]
                size = np.concatenate((size, size_tmp), axis=0)

        self.vertices = self.vertices + size

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order)


class Box_Cylinder_Ring(GeometryTemplate):
    def __init__(self, outer_height, outer_length, outer_width, inner_radius, inner_cylinder_offset = [0, 0], position = [0, 0, 0], rotation = [0, 0, 0], rotation_order = "XYZ"):
        """
        :param outer_height: height of the outer part of the box cylinder ring in the Y-axis direction
        :param outer_length: length of the outer part of the box cylinder ring in the X-axis direction
        :param outer_width: width of the outer part of the box cylinder ring in the Z-axis direction
        :param inner_radius: radius of the inner part of the box cylinder ring
        :param inner_cylinder_offset: offset between the outer and inner part of the box cylinder ring in the X-axis and Z-axis directions
        :param position: position (x, y, z) of the box cylinder ring
        :param rotation: rotation of the box cylinder ring, represented via Euler angles (x, y, z)
        :param rotation_order: rotation order of the three rotation axes of the box cylinder ring

        """

        super().__init__(position, rotation, rotation_order)

        # Record Parameters
        self.outer_height = outer_height
        self.outer_length = outer_length
        self.outer_width = outer_width
        self.inner_radius = inner_radius
        self.inner_cylinder_offset = inner_cylinder_offset
        self.position = position
        self.rotation = rotation
        self.rotation_order = rotation_order
            
        # Manually Defined Default Template Instance 
        vertices = []
        num_of_segment = 256
        vertices.append([-1 / 2, 1 / 2, 1 / 2]) # outer
        vertices.append([1 / 2, 1 / 2, 1 / 2])
        vertices.append([-1 / 2, 1 / 2, -1 / 2])
        vertices.append([1 / 2, 1 / 2, -1 / 2])
        vertices.append([-1 / 2, -1 / 2, 1 / 2])
        vertices.append([1 / 2, -1 / 2, 1 / 2])
        vertices.append([-1 / 2, -1 / 2, -1 / 2])
        vertices.append([1 / 2, -1 / 2, -1 / 2])
        
        for i in range(num_of_segment+1): # inner
            rotation_tmp = np.pi * 2 / num_of_segment * i
            vertices.append([np.cos(rotation_tmp), np.sin(rotation_tmp), 1 / 2])
            vertices.append([np.cos(rotation_tmp), np.sin(rotation_tmp), -1 / 2])
        self.vertices = np.array(vertices)

        faces = []
        faces.append([0, 1, 2]) # outer_other
        faces.append([1, 3, 2])
        faces.append([4, 6, 5])
        faces.append([5, 6, 7])
        faces.append([0, 6, 4])
        faces.append([0, 2, 6])
        faces.append([1, 5, 7])
        faces.append([1, 7, 3])

        faces.append([1, 0, 8 + num_of_segment / 4 * 1 * 2]) # outer_front
        faces.append([0, 4, 8 + num_of_segment / 4 * 2 * 2])
        faces.append([4, 5, 8 + num_of_segment / 4 * 3 * 2])
        faces.append([5, 1, 8 + num_of_segment / 4 * 4 * 2])
        label = [1, 0, 4, 5]
        for i in range(num_of_segment):
            faces.append([label[int(i * 4 / num_of_segment)], 8 + (i+1) * 2, 8 + i * 2])

        faces.append([2, 3, 8 + num_of_segment / 4 * 1 * 2 + 1]) # outer_behind
        faces.append([6, 2, 8 + num_of_segment / 4 * 2 * 2 + 1])
        faces.append([7, 6, 8 + num_of_segment / 4 * 3 * 2 + 1])
        faces.append([3, 7, 8 + num_of_segment / 4 * 4 * 2 + 1])
        label = [3, 2, 6, 7]
        for i in range(num_of_segment):
            faces.append([label[int(i * 4 / num_of_segment)], 8 + i * 2 + 1, 8 + (i+1) * 2 + 1])

        for i in range(num_of_segment): # inner
            faces.append([8 + i * 2, 8 + (i+1) * 2, 8 + i * 2 + 1])
            faces.append([8 + i * 2 + 1, 8 + (i+1) * 2, 8 + (i+1) * 2 + 1])
        self.faces = np.array(faces)

        # Differentiable Deformation
        outer_size = np.array([outer_length, outer_height, outer_width])
        outer_size = np.tile(outer_size[None, :], (8, 1))
        inner_size = np.array([inner_radius, inner_radius, outer_width])
        inner_size = np.tile(inner_size[None, :], (self.vertices.shape[0] - 8, 1))
        vertices_resize = np.concatenate((outer_size, inner_size), axis=0)
        self.vertices = self.vertices * vertices_resize

        outer_offset = np.array([0, 0, 0])
        outer_offset = np.tile(outer_offset[None, :], (8, 1))
        inner_offset = np.array([inner_cylinder_offset[0], inner_cylinder_offset[1], 0])
        inner_offset = np.tile(inner_offset[None, :], (self.vertices.shape[0] - 8, 1))
        vertices_offset = np.concatenate((outer_offset, inner_offset), axis=0)
        self.vertices = self.vertices + vertices_offset

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order)


class Cylinder_Box_Ring(GeometryTemplate):
    def __init__(self, outer_radius, outer_height, inner_length, inner_width, inner_cuboid_offset = [0, 0], position = [0, 0, 0], rotation = [0, 0, 0], rotation_order = "XYZ"):
        """
        :param outer_radius: radius of the outer part of the cylinder box ring
        :param outer_height: height of the outer part of the cylinder box ring in the Y-axis direction
        :param inner_length: length of the inner part of the cylinder box ring in the X-axis direction
        :param inner_width: width of the inner part of the cylinder box ring in the Z-axis direction
        :param inner_cuboid_offset: offset between the outer and inner part of the cylinder box ring in the X-axis and Z-axis directions
        :param position: position (x, y, z) of the cylinder box ring
        :param rotation: rotation of the cylinder box ring, represented via Euler angles (x, y, z)
        :param rotation_order: rotation order of the three rotation axes of the cylinder box ring

        """

        super().__init__(position, rotation, rotation_order)

        # Record Parameters
        self.outer_radius = outer_radius
        self.outer_height = outer_height
        self.inner_length = inner_length
        self.inner_width = inner_width
        self.inner_cuboid_offset = inner_cuboid_offset
        self.position = position
        self.rotation = rotation
        self.rotation_order = rotation_order
            
        # Manually Defined Default Template Instance 
        vertices = []
        num_of_segment = 256
        vertices.append([-1 / 2, 1 / 2, 1 / 2]) # inner
        vertices.append([1 / 2, 1 / 2, 1 / 2])
        vertices.append([-1 / 2, 1 / 2, -1 / 2])
        vertices.append([1 / 2, 1 / 2, -1 / 2])
        vertices.append([-1 / 2, -1 / 2, 1 / 2])
        vertices.append([1 / 2, -1 / 2, 1 / 2])
        vertices.append([-1 / 2, -1 / 2, -1 / 2])
        vertices.append([1 / 2, -1 / 2, -1 / 2])
        
        for i in range(num_of_segment+1): # outer
            rotation_tmp = np.pi * 2 / num_of_segment * i
            vertices.append([np.cos(rotation_tmp), np.sin(rotation_tmp), 1 / 2])
            vertices.append([np.cos(rotation_tmp), np.sin(rotation_tmp), -1 / 2])
        self.vertices = np.array(vertices)

        faces = []
        faces.append([0, 2, 1]) # outer_other
        faces.append([1, 2, 3])
        faces.append([4, 5, 6])
        faces.append([5, 7, 6])
        faces.append([0, 4, 6])
        faces.append([0, 6, 2])
        faces.append([1, 7, 5])
        faces.append([1, 3, 7])

        faces.append([0, 1, 8 + num_of_segment / 4 * 1 * 2]) # outer_front
        faces.append([4, 0, 8 + num_of_segment / 4 * 2 * 2])
        faces.append([5, 4, 8 + num_of_segment / 4 * 3 * 2])
        faces.append([1, 5, 8 + num_of_segment / 4 * 4 * 2])
        label = [1, 0, 4, 5]
        for i in range(num_of_segment):
            faces.append([label[int(i * 4 / num_of_segment)], 8 + i * 2, 8 + (i+1) * 2])

        faces.append([3, 2, 8 + num_of_segment / 4 * 1 * 2 + 1]) # outer_behind
        faces.append([2, 6, 8 + num_of_segment / 4 * 2 * 2 + 1])
        faces.append([6, 7, 8 + num_of_segment / 4 * 3 * 2 + 1])
        faces.append([7, 3, 8 + num_of_segment / 4 * 4 * 2 + 1])
        label = [3, 2, 6, 7]
        for i in range(num_of_segment):
            faces.append([label[int(i * 4 / num_of_segment)], 8 + (i+1) * 2 + 1, 8 + i * 2 + 1])

        for i in range(num_of_segment): # inner
            faces.append([8 + (i+1) * 2, 8 + i * 2, 8 + i * 2 + 1])
            faces.append([8 + (i+1) * 2 + 1, 8 + (i+1) * 2, 8 + i * 2 + 1])
        self.faces = np.array(faces)

        # Differentiable Deformation
        inner_size = np.array([inner_length, inner_width, outer_height])
        inner_size = np.tile(inner_size[None, :], (8, 1))
        outer_size = np.array([outer_radius, outer_radius, outer_height])
        outer_size = np.tile(outer_size[None, :], (self.vertices.shape[0] - 8, 1))
        vertices_resize = np.concatenate((inner_size, outer_size), axis=0)
        self.vertices = self.vertices * vertices_resize

        inner_offset = np.array([inner_cuboid_offset[0], inner_cuboid_offset[1], 0])
        inner_offset = np.tile(inner_offset[None, :], (8, 1))
        outer_offset = np.array([0, 0, 0])
        outer_offset = np.tile(outer_offset[None, :], (self.vertices.shape[0] - 8, 1))
        vertices_offset = np.concatenate((inner_offset, outer_offset), axis=0)
        self.vertices = self.vertices + vertices_offset

        # Global Transformation
        self.vertices = apply_transformation(self.vertices, position, rotation, rotation_order)