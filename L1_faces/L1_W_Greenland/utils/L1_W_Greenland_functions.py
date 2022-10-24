
import numpy as np
import matplotlib.pyplot as plt


def read_grid_faces_to_stitched_grid(face_grids, dim):


    if dim==2:
        stitched_grid = np.concatenate([np.rot90(face_grids[5], k=1),
                                        np.rot90(face_grids[3], k=2),],axis=0)
    if dim==3:
        stitched_grid = np.concatenate([np.rot90(face_grids[5],axes=(1,2), k=1),
                                        np.rot90(face_grids[3],axes=(1,2), k=2),],axis=1)


    return(stitched_grid)


def read_compact_grid_to_stitched_grid(compact_grid, sNx, sNy, faces, face_size_dict, dim):

    if dim==2:
        face_5 = compact_grid[3 * sNx:, :]
        face_5 = np.reshape(face_5, (3 * sNy, 3 * sNx))
        # plt.imshow(face_1, origin='lower')
        # plt.show()

        face_3 = compact_grid[:3*sNx,:]
        face_3 = np.reshape(face_3,(sNx,3*sNx))
        # plt.imshow(face_3, origin='lower')
        # plt.show()

        stitched_grid = np.concatenate([np.rot90(face_5, k=1),
                                        np.rot90(face_3, k=2)],axis=0)
        # stitched_grid = face_3

        # plt.imshow(stitched_grid, origin='lower')
        # plt.show()

    # if dim == 3:
    #     face_1 = compact_grid[:,:3 * sNx, :]
    #     face_1 = np.reshape(face_1, (np.shape(compact_grid)[0], sNx, 3 * sNx))
    #     # plt.imshow(face_1, origin='lower')
    #     # plt.show()
    #
    #     face_3 = compact_grid[:, 3 * sNx:, :]
    #     # plt.imshow(face_3, origin='lower')
    #     # plt.show()
    #
    #     stitched_grid = np.concatenate([face_1, np.rot90(face_3, k=3, axes=(1,2))],axis=1)
    #     # plt.imshow(stitched_grid, origin='lower')
    #     # plt.show()

    if dim==3:
        face_5 = compact_grid[:, 3 * sNx:, :]
        face_5 = np.reshape(face_5, (np.shape(compact_grid)[0], 3 * sNy, 3 * sNx))
        # plt.imshow(face_1, origin='lower')
        # plt.show()

        face_3 = compact_grid[:,:3*sNx,:]
        face_3 = np.reshape(face_3,(np.shape(compact_grid)[0], sNx,3*sNx))
        # plt.imshow(face_3, origin='lower')
        # plt.show()

        stitched_grid = np.concatenate([np.rot90(face_5, k=1, axes=(1,2)),
                                        np.rot90(face_3, k=2, axes=(1,2))],axis=1)


    return(stitched_grid)


def read_stitched_grid_to_faces(stitched_grid, sNx, sNy, dim):

    face_grids = {}

    if dim == 2:
        face_grids[3] = np.rot90(stitched_grid[-sNy:,:],k=2)
        face_grids[5] = np.rot90(stitched_grid[:3*sNy,:],k=3)

    if dim == 3:
        face_grids[3] = np.rot90(stitched_grid[:,-sNy:,:],axes = (1,2),k=2)
        face_grids[5] = np.rot90(stitched_grid[:,:3*sNy,:],axes = (1,2),k=3)

    return(face_grids)


def read_faces_to_compact(grid_faces, sNx, sNy):

    faces = list(grid_faces.keys())

    for f in range(len(faces)):
        face = faces[f]
        grid = grid_faces[face]
        n_rows = int((np.shape(grid)[1] * np.shape(grid)[2]) / sNx)
        grid = np.reshape(grid, (np.shape(grid)[0], n_rows, sNx))
        if f==0:
            compact_stack = grid
        else:
            compact_stack = np.concatenate([compact_stack, grid], axis=1)

    return(compact_stack)

