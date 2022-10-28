
import numpy as np
import matplotlib.pyplot as plt


def read_grid_faces_to_stitched_grid(face_grids, dim):


    if dim==2:
        stitched_grid = np.concatenate([np.rot90(face_grids[5], k=1),
                                        face_grids[1]],axis=1)
    if dim==3:
        stitched_grid = np.concatenate([np.rot90(face_grids[5],axes=(1,2), k=1),
                                        face_grids[1]],axis=2)


    return(stitched_grid)


def read_compact_grid_to_stitched_grid(compact_grid, sNx, sNy, faces, face_size_dict, dim):

    if dim==2:
        face_5 = compact_grid[-2 * sNx:, :]
        face_5 = np.reshape(face_5, (sNy, 2 * sNx))
        # plt.imshow(face_1, origin='lower')
        # plt.show()

        face_1 = compact_grid[:4*sNx,:]
        face_1 = np.reshape(face_1,(2*sNy,2*sNx))
        # plt.imshow(face_3, origin='lower')
        # plt.show()

        stitched_grid = np.concatenate([np.rot90(face_5, k=1),
                                        face_1],axis=1)
        # stitched_grid = face_3

        # plt.imshow(stitched_grid, origin='lower')
        # plt.show()

    if dim == 3:
        face_5 = compact_grid[:, -2 * sNx:, :]
        face_5 = np.reshape(face_5, (np.shape(compact_grid)[0], sNy, 2 * sNx))

        face_1 = compact_grid[:,:4 * sNx, :]
        face_1 = np.reshape(face_1, (np.shape(compact_grid)[0], 2* sNx, 2 * sNx))
        # plt.imshow(face_1, origin='lower')
        # plt.show()

        stitched_grid = np.concatenate([np.rot90(face_5, k=1, axes=(1,2)), face_1],axis=2)
        # plt.imshow(stitched_grid, origin='lower')
        # plt.show()


    return(stitched_grid)


def read_stitched_grid_to_faces(stitched_grid, sNx, sNy, dim):

    face_grids = {}

    if dim == 2:
        face_grids[1] = stitched_grid[:,sNx:]
        face_grids[5] = np.rot90(stitched_grid[:, :sNx], k=3)

    if dim == 3:
        face_grids[1] = stitched_grid[:,sNx:,:]
        face_grids[5] = np.rot90(stitched_grid[:,:sNx,:],axes = (1,2), k=3)

    return(face_grids)


def read_faces_to_compact(grid_faces, sNx, sNy):

    for face in list(grid_faces.keys()):
        grid = grid_faces[face]
        n_rows = int((np.shape(grid)[1] * np.shape(grid)[2]) / sNx)
        grid = np.reshape(grid, (np.shape(grid)[0], n_rows, sNx))
        if face == 1:
            compact_stack = grid
        else:
            compact_stack = np.concatenate([compact_stack, grid], axis=1)

    return(compact_stack)



