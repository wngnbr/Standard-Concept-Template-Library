from concept_template import *
from geometry_template import *
from scipy.spatial.transform import Rotation as Rot
import pickle
from utils import COLOR20
import open3d as o3d
from knowledge_utils import *
from knowledge_definitions import *


if __name__ == "__main__":

    with open("conceptualization.pkl", "rb") as f:
        data_list = pickle.load(f)

    for data in data_list:

        conceptualization = data['conceptualization']

        objs = {}
        region_knowledge_wrappers = {}


        vertices_list = []
        faces_list = []
        total_num_vertices = 0
        
        for c in conceptualization:
            module = eval(c['template'])
            obj = module(**c['parameters'])
            vertices_list.append(obj.vertices)
            faces_list.append(obj.faces + total_num_vertices)
            total_num_vertices += len(obj.vertices)

            objs[c['template']] = obj
            region_knowledge_wrappers[c['template']] = Region_Knowledge_Wrapper(obj)

        final_vertices = np.concatenate(vertices_list)
        final_faces = np.concatenate(faces_list)
        final_mesh = trimesh.Trimesh(final_vertices, final_faces)
        pts = np.array(final_mesh.sample(20000))
        

        # affordance
        affordance_label = np.zeros(pts.shape[0])
        for template, obj in objs.items():
            res = np.full((pts.shape[0]), False)
            if (obj.semantic == 'Lid'):
                res = region_knowledge_wrappers[template].check(lid_affordance, pts)
            for idx in range(pts.shape[0]):
                if (affordance_label[idx] == 0 and res[idx] != False):
                    affordance_label[idx] = res[idx]

        affordance_label = affordance_label.astype(np.int32)
        color = COLOR20[affordance_label]
        pcd = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(pts))
        pcd.colors = o3d.utility.Vector3dVector(color)
        o3d.visualization.draw_geometries([pcd])


        # part pose
        coordinates = []
        poses = []
        for template, obj in objs.items():

            pose = part_pose(obj)
            poses.append(poses)

            coordinate = o3d.geometry.TriangleMesh.create_coordinate_frame(size=0.5)
            v = np.array(coordinate.vertices)
            v = np.concatenate([v, np.ones((v.shape[0], 1))], axis=1)
            v = (pose @ v.T).T
            coordinate.vertices = o3d.utility.Vector3dVector(v[:, :3])
            coordinates.append(coordinate)
        
        pcd = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(pts))
        o3d.visualization.draw_geometries([pcd] + coordinates)