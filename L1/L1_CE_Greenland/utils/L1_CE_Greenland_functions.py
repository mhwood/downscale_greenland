
import numpy as np
import matplotlib.pyplot as plt


def read_grid_faces_to_stitched_grid(face_grids, dim):

    if dim==3:
        # stitched_grid = np.concatenate([face_grids[1], np.transpose(face_grids[3],axes=(0,2,1))],axis=1)
        stitched_grid = np.concatenate([face_grids[1],
                                        np.flip(np.flip(np.rot90(face_grids[3], axes=(1,2)),axis=2), axis=1)], axis=1)
        # stitched_grid = np.concatenate([np.flip(np.rot90(face_grids[3], axes=(1, 2)), axis=2),
        #                                 face_grids[1]], axis=1)

    return(stitched_grid)


def read_compact_grid_to_stitched_grid(compact_grid, sNx, sNy, faces, face_size_dict):

    face_1 = compact_grid[:3 * sNx, :]
    face_1 = np.reshape(face_1, (sNx, 3 * sNx))
    # plt.imshow(face_1, origin='lower')
    # plt.show()

    face_3 = compact_grid[3*sNx:,:]
    # plt.imshow(face_3, origin='lower')
    # plt.show()

    stitched_grid = np.concatenate([face_1,np.rot90(face_3,k=3)])
    # plt.imshow(stitched_grid, origin='lower')
    # plt.show()

    return(stitched_grid)








