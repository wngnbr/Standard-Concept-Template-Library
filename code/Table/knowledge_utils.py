import numpy as np
from scipy.spatial.transform import Rotation as Rot


PROXIMITY_THRES = 0.02
AFFORDACE_PROXIMITY_THRES = 0.02
SAMPLENUM = 10000


class Region_Knowledge_Wrapper():
    def __init__(self, template_instance):
        self.instance = template_instance
    
    def check(self, func, pts, *args, **kwargs):
        res = []
        for pt in pts:
            if (not self.instance.proximation(pt)):
                res.append(False)
            else:
                res.append(func(self.instance, pt, *args, **kwargs))

        return res


def transformation_matrix(position, rotation):
    RT = np.eye(4)
    RT[:3, :3] = Rot.from_euler('xyz', rotation, degrees=False).as_matrix()
    RT[:3, -1] = position
    return RT


def inverse_transformation(pt, position, rotation):
    RT = transformation_matrix(position, rotation)
    _pt = np.array([pt[0], pt[1], pt[2], 1])
    _pt = np.linalg.inv(RT) @ _pt
    return _pt[:3]