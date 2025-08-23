#!/usr/bin/env python

import numpy as np

"""
Based on a test (see main code below), using quaternion is about twice faster
than Avro's method. Therefore, using quaternion is the default.
"""

"""
Algorithm by Avro 1992
"""
def uniform_random_rotation_Avro(x):
    # (c) Jack Scantlebury
    # https://www.blopig.com/blog/2021/08/uniformly-sampled-3d-rotation-matrices/
    """Apply a random rotation in 3D, with a distribution uniform over the
    sphere.

    Arguments:
        x: vector or set of vectors with dimension (n, 3), where n is the
            number of vectors

    Returns:
        Array of shape (n, 3) containing the randomly rotated vectors of x,
        about the mean coordinate of x.

    Algorithm taken from "Fast Random Rotation Matrices" (James Avro, 1992):
    https://doi.org/10.1016/B978-0-08-050755-2.50034-8
    """

    def generate_random_z_axis_rotation():
        """Generate random rotation matrix about the z axis."""
        R = np.eye(3)
        x1 = np.random.rand()
        R[0, 0] = R[1, 1] = np.cos(2 * np.pi * x1)
        R[0, 1] = -np.sin(2 * np.pi * x1)
        R[1, 0] = np.sin(2 * np.pi * x1)
        return R

    # There are two random variables in [0, 1) here (naming is same as paper)
    x2 = 2 * np.pi * np.random.rand()
    x3 = np.random.rand()

    # Rotation of all points around x axis using matrix
    R = generate_random_z_axis_rotation()
    v = np.array([
        np.cos(x2) * np.sqrt(x3),
        np.sin(x2) * np.sqrt(x3),
        np.sqrt(1 - x3)
    ])
    H = np.eye(3) - (2 * np.outer(v, v))
    M = -(H @ R)
    x = np.array(x).reshape((-1, 3))
    mean_coord = np.mean(x, axis=0)

    return ((x - mean_coord) @ M) + mean_coord @ M

"""
Algorithm using random quaternion
"""
def random_rotation_matrix_by_quaternion():
    """
    Generate a uniformly random rotation matrix
    """
    u1, u2, u3 = np.random.rand(3)
    q1 = np.sqrt(1 - u1) * np.sin(2 * np.pi * u2)
    q2 = np.sqrt(1 - u1) * np.cos(2 * np.pi * u2)
    q3 = np.sqrt(u1) * np.sin(2 * np.pi * u3)
    q4 = np.sqrt(u1) * np.cos(2 * np.pi * u3)

    # From quaternion to a rotation matrix
    r11 = 1 - 2*(q3**2 + q4**2)
    r12 = 2*(q2*q3 - q1*q4)
    r13 = 2*(q2*q4 + q1*q3)
    r21 = 2*(q2*q3 + q1*q4)
    r22 = 1 - 2*(q2**2 + q4**2)
    r23 = 2*(q3*q4 - q1*q2)
    r31 = 2*(q2*q4 - q1*q3)
    r32 = 2*(q3*q4 + q1*q2)
    r33 = 1 - 2*(q2**2 + q3**2)
    return np.array([[r11, r12, r13],
                     [r21, r22, r23],
                     [r31, r32, r33]])

def uniform_random_rotation(x):
    R = random_rotation_matrix_by_quaternion()
    x = np.array(x).reshape((-1, 3))
    return np.dot(x, R.T)


if __name__ == "__main__":

    import time

    N = 1000000
    x = [[0,0,0],[10,0,0],[0,10,0]]
    ## Generate random rotation of this object x for N times.
    # Output the coordinates in XYZ format.
    # Random distribution can be visually confirmed by plotting the points in 3D,
    # or analyse the coordinates.
    # Also measure the time taken by the two methods.

    """ Avro's method """
    #fout0 = open('test_Avro_0.xyz', 'w')  # This point is always [0, 0, 0]
    #fout1 = open('test_Avro_1.xyz', 'w')
    fout2 = open('test_Avro_2.xyz', 'w')
    #fout0.write(f'{N}\n\n')
    #fout1.write(f'{N}\n\n')
    fout2.write(f'{N}\n\n')

    start = time.time()

    for i in range(N):
        y = uniform_random_rotation_Avro(x)
        #fout0.write('C %f %f %f\n' % (y[0][0], y[0][1], y[0][2]))
        #fout1.write('C %f %f %f\n' % (y[1][0], y[1][1], y[1][2]))
        fout2.write('C %f %f %f\n' % (y[2][0], y[2][1], y[2][2]))

    end = time.time()
    print('Avro 1992: ', end - start)

    #fout0.close()
    #fout1.close()
    fout2.close()


    """ Using quaternion """
    #foutQ0 = open('test_quaternion_0.xyz', 'w')
    #foutQ1 = open('test_quaternion_1.xyz', 'w')
    foutQ2 = open('test_quaternion_2.xyz', 'w')

    #foutQ0.write(f'{N}\n\n')
    #foutQ1.write(f'{N}\n\n')
    foutQ2.write(f'{N}\n\n')

    start = time.time()
    for i in range(N):
        y = uniform_random_rotation(x)
        #foutQ0.write('C %f %f %f\n' % (y[0][0], y[0][1], y[0][2]))
        #foutQ1.write('C %f %f %f\n' % (y[1][0], y[1][1], y[1][2]))
        foutQ2.write('C %f %f %f\n' % (y[2][0], y[2][1], y[2][2]))
    end = time.time()
    print('Quaternion: ', end - start)

    #foutQ0.close()
    #foutQ1.close()
    foutQ2.close()
