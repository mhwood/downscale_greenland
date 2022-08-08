
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


def read_stitched_grid_to_faces(stitched_grid, sNx, sNy, dim):

    face_grids = {}

    if dim == 2:
        face_grids[1] = stitched_grid[:sNy,:]
        face_grids[3] = np.rot90(stitched_grid[sNy:,:])

    if dim == 3:
        face_grids[1] = stitched_grid[:,:sNy,:]
        face_grids[3] = np.rot90(stitched_grid[:,sNy:,:],axes = (1,2))

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


def read_compact_grid_to_stitched_grid(compact_grid, sNx, sNy, faces, face_size_dict, dim):

    if dim==2:
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

    if dim == 3:
        face_1 = compact_grid[:,:3 * sNx, :]
        face_1 = np.reshape(face_1, (np.shape(compact_grid)[0], sNx, 3 * sNx))
        # plt.imshow(face_1, origin='lower')
        # plt.show()

        face_3 = compact_grid[:, 3 * sNx:, :]
        # plt.imshow(face_3, origin='lower')
        # plt.show()

        stitched_grid = np.concatenate([face_1, np.rot90(face_3, k=3, axes=(1,2))],axis=1)
        # plt.imshow(stitched_grid, origin='lower')
        # plt.show()


    return(stitched_grid)










