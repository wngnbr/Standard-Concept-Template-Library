from concept_template import *
from geometry_template import *
from knowledge_utils import *


def handle_affordance(obj, pt):

    def is_affordance(obj, pt):
        if (isinstance(obj, Cuboidal_Handle)):
            affordable_handle = np.array([[-obj.size[0], -obj.size[1], 0],
                                          [obj.size[0], obj.size[1], obj.size[2]]])
            affordable_handle[0] -= AFFORDACE_PROXIMITY_THRES
            affordable_handle[1] += AFFORDACE_PROXIMITY_THRES
            
            _pt = inverse_transformation(pt, obj.position, obj.rotation)

            return (_pt[0] >= affordable_handle[0, 0] and 
                    _pt[0] <= affordable_handle[1, 0] and 
                    _pt[1] >= affordable_handle[0, 1] and 
                    _pt[1] <= affordable_handle[1, 1] and 
                    _pt[2] >= affordable_handle[0, 2] and 
                    _pt[2] <= affordable_handle[1, 2])
        
        elif (isinstance(obj, Trifold_Handle)):
            affordable_handle = np.array([[-obj.grip_size[0], -obj.grip_size[1], obj.grip_size[2]],
                                          [obj.grip_size[0], obj.grip_size[1], obj.grip_size[2] + obj.mounting_size[2]]])
            affordable_handle[0] -= AFFORDACE_PROXIMITY_THRES
            affordable_handle[1] += AFFORDACE_PROXIMITY_THRES
            
            _pt = inverse_transformation(pt, obj.position, obj.rotation)

            return (_pt[0] >= affordable_handle[0, 0] and 
                    _pt[0] <= affordable_handle[1, 0] and 
                    _pt[1] >= affordable_handle[0, 1] and 
                    _pt[1] <= affordable_handle[1, 1] and 
                    _pt[2] >= affordable_handle[0, 2] and 
                    _pt[2] <= affordable_handle[1, 2])
    
        elif (isinstance(obj, Curve_Handle)):
            _pt = inverse_transformation(pt, obj.position, obj.rotation)

            curve_z_offset = -obj.curve_size[0] * np.cos(obj.curve_exist_angle[0] / 2)
            _pt_radius = np.linalg.norm(_pt - np.array([0, 0, -obj.curve_size[0] * np.cos(obj.curve_exist_angle[0] / 2)]))
            _pt_angle = np.arccos((_pt[2]-curve_z_offset) / np.sqrt(_pt[1]**2 + (_pt[2]-curve_z_offset)**2))

            return (_pt_radius >= obj.curve_size[1] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.curve_size[0] + AFFORDACE_PROXIMITY_THRES and
                    np.abs(_pt[0]) <= obj.curve_size[2]/2 and 
                    _pt_angle >= 0 and _pt_angle <= obj.curve_exist_angle[0]/2)
        

        elif (isinstance(obj, Trifold_Curve_Handle)):
            _pt = inverse_transformation(pt, obj.position, obj.rotation)

            curve_z_offset = obj.mounting_size[2] - np.sqrt(obj.curve_size[1]**2 - (obj.mounting_seperation[0] / 2)**2)
            _pt_radius = np.linalg.norm(_pt - np.array([0, 0, curve_z_offset]))
            _pt_angle = np.arccos((_pt[2]-curve_z_offset) / np.sqrt(_pt[1]**2 + (_pt[2]-curve_z_offset)**2))
            
            return (_pt_radius >= obj.curve_size[1] - AFFORDACE_PROXIMITY_THRES and 
                    _pt_radius <= obj.curve_size[0] + AFFORDACE_PROXIMITY_THRES and
                    np.abs(_pt[0]) <= obj.curve_size[2]/2 and 
                    _pt_angle >= 0 and _pt_angle <= obj.curve_exist_angle[0]/2)
    
    
    return is_affordance(obj, pt)



def cuboidal_door_center(obj, pt):
    door_init_bounding_box = np.array([[-obj.size[0]/2, -obj.size[1]/2, 0],
                                        [obj.size[0]/2, obj.size[1]/2, obj.size[2]]])
    door_init_bounding_box[0] -= PROXIMITY_THRES
    door_init_bounding_box[1] += PROXIMITY_THRES

    _pt = inverse_transformation(pt, obj.position, obj.rotation)

    return (_pt[0] >= door_init_bounding_box[0, 0]*2/3 + door_init_bounding_box[1, 0]*1/3 and 
            _pt[0] <= door_init_bounding_box[0, 0]*1/3 + door_init_bounding_box[1, 0]*2/3 and 
            _pt[1] >= door_init_bounding_box[0, 1]*2/3 + door_init_bounding_box[1, 1]*1/3 and 
            _pt[1] <= door_init_bounding_box[0, 1]*1/3 + door_init_bounding_box[1, 1]*2/3 and 
            _pt[2] >= door_init_bounding_box[0, 2]*2/3 + door_init_bounding_box[1, 2]*1/3 and 
            _pt[2] <= door_init_bounding_box[0, 2]*1/3 + door_init_bounding_box[1, 2]*2/3)


def sunken_door_center(obj, pt):
    door_init_bounding_box = np.array([[-obj.size[0]/2, -obj.size[1]/2, 0],
                                       [obj.size[0]/2, obj.size[1]/2, obj.size[2]-obj.sunken_size[2]]])
    door_init_bounding_box[0] -= PROXIMITY_THRES
    door_init_bounding_box[1] += PROXIMITY_THRES

    _pt = inverse_transformation(pt, obj.position, obj.rotation)

    return (_pt[0] >= door_init_bounding_box[0, 0]*2/3 + door_init_bounding_box[1, 0]*1/3 and 
            _pt[0] <= door_init_bounding_box[0, 0]*1/3 + door_init_bounding_box[1, 0]*2/3 and 
            _pt[1] >= door_init_bounding_box[0, 1]*2/3 + door_init_bounding_box[1, 1]*1/3 and 
            _pt[1] <= door_init_bounding_box[0, 1]*1/3 + door_init_bounding_box[1, 1]*2/3 and 
            _pt[2] >= door_init_bounding_box[0, 2]*2/3 + door_init_bounding_box[1, 2]*1/3 and 
            _pt[2] <= door_init_bounding_box[0, 2]*1/3 + door_init_bounding_box[1, 2]*2/3)



def part_pose(obj):
    RT = transformation_matrix(obj.position, obj.rotation)
    return RT



def grasp_poses(obj, num):
    grasp_poses = []

    if (isinstance(obj, Cuboidal_Handle)):
        RT = transformation_matrix(obj.position, obj.rotation)
        
        for height_ratio in np.random.random(num):
            init_pose = np.eye(4)
            init_pose[:3, -1] = np.array([0,
                                          (height_ratio - 0.5) * obj.size[1],
                                          obj.size[2]])
            grasp_pose = RT @ init_pose
            grasp_poses.append(grasp_pose)

    elif (isinstance(obj, Trifold_Handle)):
        RT = transformation_matrix(obj.position, obj.rotation)
        
        for height_ratio in np.random.random(num):
            init_pose = np.eye(4)
            init_pose[:3, -1] = np.array([0,
                                          (height_ratio - 0.5) * obj.top_size[1],
                                          obj.bottom_size[2] + obj.top_size[2]])
            grasp_pose = RT @ init_pose
            grasp_poses.append(grasp_pose)

    elif (isinstance(obj, Curve_Handle)):
        RT = transformation_matrix(obj.position, obj.rotation)

        curve_z_offset = obj.bottom_size[2] - np.sqrt(obj.curve_size[1]**2 - (obj.bottom_seperation[0] / 2)**2)
        
        for angle_ratio in np.random.random(num):
            init_pose = np.eye(4)
            init_pose[:3, -1] = np.array([0,
                                          np.sin(angle_ratio - 0.5) * obj.curve_size[0],
                                          np.cos(angle_ratio - 0.5) * obj.curve_size[0] + curve_z_offset])
            grasp_pose = RT @ init_pose
            grasp_poses.append(grasp_pose)

    elif (isinstance(obj, Trifold_Curve_Handle)):
        RT = transformation_matrix(obj.position, obj.rotation)
        
        curve_z_offset = obj.bottom_size[2] - np.sqrt(obj.curve_size[1]**2 - (obj.bottom_seperation[0] / 2)**2)

        for angle_ratio in np.random.random(num):
            init_pose = np.eye(4)
            init_pose[:3, -1] = np.array([0,
                                          np.sin(angle_ratio - 0.5) * obj.curve_size[0],
                                          np.cos(angle_ratio - 0.5) * obj.curve_size[0] + curve_z_offset])
            grasp_pose = RT @ init_pose
            grasp_poses.append(grasp_pose)
    
    return grasp_poses