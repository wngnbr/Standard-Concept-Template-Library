import numpy as np
import os
import open3d as o3d



COLOR20 = np.array(
    [[230, 230, 230], [0, 128, 128], [230, 190, 255], [170, 110, 40], [255, 250, 200], [128, 0, 0],
     [170, 255, 195], [128, 128, 0], [255, 215, 180], [0, 0, 128], [128, 128, 128],
     [230, 25, 75], [60, 180, 75], [255, 225, 25], [0, 130, 200], [245, 130, 48],
     [145, 30, 180], [70, 240, 240], [240, 50, 230], [210, 245, 60], [250, 190, 190]]) / 255


def get_rodrigues_matrix(axis, angle):
    axis = np.array(axis)
    identity = np.eye(3)
    s1 = np.array(
        [
            [np.zeros([]), -axis[2], axis[1]],
            [axis[2], np.zeros([]), -axis[0]],
            [-axis[1], axis[0], np.zeros([])],
        ]
    )
    s2 = np.matmul(axis[:, None], axis[None])
    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)

    rodrigues_matrix = cos_angle * identity + sin_angle * s1 + (1 - cos_angle) * s2
    return rodrigues_matrix


def apply_transformation(vertices, position, rotation, rotation_order="XYZ", offset_first=False):

    # process position first
    if offset_first:
        vertices = vertices + np.array(position)

    # process rotation
    rot_mat = {}

    rot_mat["X"] = get_rodrigues_matrix([1, 0, 0], rotation[0])
    rot_mat["Y"] = get_rodrigues_matrix([0, 1, 0], rotation[1])
    rot_mat["Z"] = get_rodrigues_matrix([0, 0, 1], rotation[2])

    for s in rotation_order:
        vertices = np.matmul(vertices, rot_mat[s].T)

    # process position second
    if not offset_first:
        vertices = vertices + np.array(position)

    return vertices


def adjust_position_from_rotation(position, rotation, rotation_order="XYZ"):

    position_new = np.array(position)[None, :]
    position_new = apply_transformation(
        position_new, [0, 0, 0], rotation, rotation_order
    )

    return list(position_new[0])


def list_add(list1, list2):
    list_res = []
    for i in range(len(list1)):
        list_res.append(list1[i] + list2[i])
    return list_res
