import numpy as np

def numerical(theta_2, theta_3, theta_t, num_points=int(1e6)):
    # This is a simplified version of the script written by
    # Andrew G. York "spherical_cap_overlap.py" from:
    # https://andrewgyork.github.io/high_na_single_objective_lightsheet/
    # appendix.html#theoretical_performance
    
    # 2 spherical caps: 2 centered on the z-axis, 3 tilts over by theta_t
    
    # Generate random points on the surface of a sphere:
    phi = np.random.uniform(0, 2 * np.pi, num_points)
    # https://en.wikipedia.org/wiki/Inverse_transform_sampling
    theta = np.arccos(np.random.uniform(-1, 1, num_points))
    
    # Calulate z in a rotated frame where cap 2 is centered on the (new) z-axis:
    # Dot product of unit vectors to transform from z to z_rotated
    # may help: https://en.wikipedia.org/wiki/Rotation_matrix
    z_rotated = (np.sin(theta_t) * np.sin(theta) * np.sin(phi) +
                 np.cos(theta_t) * np.cos(theta))
    
    # Check which caps each point occupies:
    in_cap2 = theta < theta_2
    in_cap3 = z_rotated > np.cos(theta_3)
    
    # Estimate the fraction of the 2nd cap covered by the 3rd cap
    ratio_2_3 = np.count_nonzero(in_cap2 & in_cap3) / np.count_nonzero(in_cap2)
    return ratio_2_3

def analytical(theta_2, theta_3, theta_t):
    if theta_t >= theta_2 + theta_3: # NA3 outside NA2 -> no transmission
        trans_3 = 0
    elif theta_t + theta_2 <= theta_3: # NA2 inside NA3 -> full transmission
        trans_3 = 1
    elif theta_t + theta_3 <= theta_2: # NA3 inside NA2 -> spherical cap ratio
        trans_3 = (1 - np.cos(theta_3)) / (1 - np.cos(theta_2))
    else: # NA3 and NA2 intersect
        # area_2_3 = overlap of spherical caps 2 and 3 seperated by theta_t:
        # "https://math.stackexchange.com/questions/45788/calculate-the-area
        # -on-a-sphere-of-the-intersection-of-two-spherical-caps"
        a, b, c, r = theta_2, theta_3, theta_t, 1 # unit radius
        s = 0.5 * (a + b + c)
        k = ((np.sin(s - a) * np.sin(s - b) * np.sin(s - c)) / np.sin(s))**0.5
        A = 2 * np.arctan(k / np.sin(s - a))
        B = 2 * np.arctan(k / np.sin(s - b))
        C = 2 * np.arctan(k / np.sin(s - c))
        area_2_3 = 2 * r**2 * (np.pi - A * np.cos(b) - B * np.cos(a) - C)
        area_2 = 2 * np.pi * r**2 * (1 - np.cos(a)) # area of spherical cap 2
        trans_3 = area_2_3 / area_2
    return trans_3

if __name__ == "__main__":
    print('Intersection of spherical caps:')

    # single value:
    theta_2, theta_3, theta_t = np.deg2rad((72, 90, 30)) # ~ SOLS config
    ratio_2_3 = numerical(theta_2, theta_3, theta_t)
    trans_3 = analytical(theta_2, theta_3, theta_t)
    print('numerical  = %0.2f%%'%(100 * ratio_2_3))
    print('analytical = %0.2f%%'%(100 * trans_3))
    error = 100 * (trans_3 - ratio_2_3)
    print('error = %0.2f%%\n'%error)

    # error check:
    num_points = int(1e5) # ~stable at 1e6 for 100 iterations...
    iterations = 10
    for i in range(iterations):
        theta_2, theta_3, theta_t = np.random.uniform(0, np.pi / 2, 3)
        trans_3 = analytical(theta_2, theta_3, theta_t)
        ratio_2_3 = numerical(theta_2, theta_3, theta_t, num_points=num_points)    
        error = 100 * (trans_3 - ratio_2_3)
        if error > 1:
            print('error > 1%% !!!')
            raise
        else:
            print('passed iteration: %2d'%i)
    print('error < 1%% (num_points=%i, iterations=%i)'%(num_points, iterations))
    
