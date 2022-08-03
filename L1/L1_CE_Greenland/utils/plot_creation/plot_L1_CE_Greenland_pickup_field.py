import os
import argparse
import ast
import sys


def plot_pickup(config_dir):
    sys.path.insert(1, os.path.join(config_dir, 'L1', 'utils', 'plot_creation'))
    import plot_pickup_field_faces as pf

    model_name = 'L1_CE_Greenland'
    llc = 1080

    sNx = 180
    sNy = 180

    tile_face_index_dict = {1: [1, 0, 0],
                            2: [1, 0, sNx],
                            3: [1, 0, 2 * sNx],
                            4: [3, 0, 0],
                            5: [3, sNy, 0],
                            6: [3, 2 * sNy, 0]}

    ordered_nonblank_tiles = [[1, 2, 3], [6, 5, 4]]

    faces_dim_file = os.path.join(config_dir, 'L1', model_name, 'namelist', 'face_dimensions.txt')
    f = open(faces_dim_file)
    dict_str = f.read()
    f.close()
    size_dict = ast.literal_eval(dict_str)
    faces = list(size_dict.keys())

    field_name = 'Theta'

    pf.plot_pickup_field(config_dir, model_name, field_name, faces, size_dict)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_pickup(config_dir)


