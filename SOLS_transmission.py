import numpy as np

def obj1_ray_transmission(theta_1):
    trans_1 = 1 - np.cos(theta_1) # fraction of hemisphere
    return trans_1

def obj2_ray_transmission(theta_1, theta_2):
    if theta_2 >= theta_1: # all rays pass
        trans_2 = 1
    else: # fraction of spherical cap 2 to spherical cap 1
        trans_2 = (1 - np.cos(theta_2)) / (1 - np.cos(theta_1))
    return trans_2

def obj3_ray_transmission(theta_2, theta_3, theta_t):
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

def ray_transmission(NA1, n1, NA2, n2, NA3, n3, theta_t_deg, verbose=True):
    theta_1 = np.arcsin(NA1 / n1)
    trans_1 = obj1_ray_transmission(theta_1)
    if verbose:
        print('objective_1: NA1=%0.2f, n1=%0.2f, '%(NA1, n1) +
              'theta_1=%0.2fdeg, trans_1=%0.2f%%'%(
                  np.rad2deg(theta_1), (100 * trans_1)))
    if NA2 is None:
        trans_2 = 1
    else:
        theta_2 = np.arcsin(NA2 / n2)
        trans_2 = obj2_ray_transmission(theta_1, theta_2)
        if verbose:
            print('objective_2: NA2=%0.2f, n2=%0.2f, '%(NA2, n2) +
                  'theta_2=%0.2fdeg, trans_2=%0.2f%%'%(
                      np.rad2deg(theta_2), (100 * trans_2)))
    if NA3 is None:
        trans_3 = 1
    else:
        if NA3 >= n2:
            theta_3 = np.pi / 2 # hemispheric collection
        else:
            theta_3 = np.arcsin(NA3 / n2) # angle of O3 in the O2 immersion
        theta_t = np.deg2rad(theta_t_deg)
        trans_3 = obj3_ray_transmission(theta_2, theta_3, theta_t)
        if verbose:
            print('objective_3: NA3=%0.2f, n3=%0.2f, '%(NA3, n3) +
                  'theta_3=%0.2fdeg, trans_3=%0.2f%%'%(
                      np.rad2deg(theta_3), (100 * trans_3)))
    ray_trans = trans_1 * trans_2 * trans_3
    if verbose:
        print('Ray transmission = %0.2f%%'%(100 * ray_trans))
    return ray_trans

def optic_transmission(
    objectives, scan_lenses, mirrors, dichroics, filters, verbose=True):
    tube_lenses = objectives # every infinity objective has a tube lens
    # approximate efficiencies hard coded for ~532nm
    opt_trans = (0.875**objectives *
                 0.975**tube_lenses *
                 0.90**scan_lenses *
                 0.99**mirrors *
                 0.95**dichroics *
                 0.95**filters)
    if verbose:
        print('Optics transmission = %0.2f%%'%(100 * opt_trans))
    return opt_trans

def total_transmission(
    name,               # system name for printed output
    objectives,         # number of objectives
    scan_lenses,        # number of scan lenses
    mirrors,            # number of mirrors
    dichroics,          # number of dichroics
    filters,            # number of emission filters
    NA1,                # objective 1 numerical aperture
    n1,                 # objective 1 immersion
    NA2,                # objective 2 numerical aperture
    n2,                 # objective 2 immersion
    NA3,                # objective 3 numerical aperture
    n3,                 # objective 3 immersion
    theta_t_deg,        # objective 2 to objective 3 tilt angle
    verbose=True):      # print output
    print('System: %s'%name)
    ray_trans = ray_transmission(
        NA1, n1, NA2, n2, NA3, n3, theta_t_deg, verbose=verbose)
    opt_trans = optic_transmission(
        objectives, scan_lenses, mirrors, dichroics, filters, verbose=verbose)
    tot_trans = ray_trans * opt_trans
    print('Total transmission = %0.2f%%'%(100 * tot_trans))
    return tot_trans

if __name__ == "__main__":
    verbose = True

    print('***Benchmark microscope:***\n')

    tot_trans = total_transmission(name='High NA epifluorescence',
                                   objectives=1,
                                   scan_lenses=0,
                                   mirrors=1,
                                   dichroics=1,
                                   filters=1,
                                   NA1=1.35,
                                   n1=1.41,
                                   NA2=None,
                                   n2=None,
                                   NA3=None,
                                   n3=None,
                                   theta_t_deg=None,
                                   verbose=verbose)
    max_trans = tot_trans
    print('Normalised transmission = %0.2f%%\n'%(100 * tot_trans / max_trans))

    print('***Dipping microscopes:***\n')

    tot_trans = total_transmission(name='SCAPE 2015, high NA config (tilt=53?)',
                                   objectives=3,
                                   scan_lenses=2,
                                   mirrors=2,
                                   dichroics=1,
                                   filters=1,
                                   NA1=1.00,
                                   n1=1.33,
                                   NA2=0.75,
                                   n2=1.00,
                                   NA3=0.40,
                                   n3=1.00,
                                   theta_t_deg=53,
                                   verbose=verbose)
    print('Normalised transmission = %0.2f%%\n'%(100 * tot_trans / max_trans))

    tot_trans = total_transmission(name='SOPi 2018',
                                   objectives=3,
                                   scan_lenses=2,
                                   mirrors=1,
                                   dichroics=2,
                                   filters=1,
                                   NA1=1.00,
                                   n1=1.33,
                                   NA2=0.75,
                                   n2=1.00,
                                   NA3=0.45,
                                   n3=1.00,
                                   theta_t_deg=45,
                                   verbose=verbose)
    print('Normalised transmission = %0.2f%%\n'%(100 * tot_trans / max_trans))

    tot_trans = total_transmission(name='SCAPE 2019, high NA config',
                                   objectives=3,
                                   scan_lenses=2,
                                   mirrors=2,
                                   dichroics=1,
                                   filters=1,
                                   NA1=1.00,
                                   n1=1.33,
                                   NA2=0.75,
                                   n2=1.00,
                                   NA3=0.75,
                                   n3=1.00,
                                   theta_t_deg=53,
                                   verbose=verbose)
    print('Normalised transmission = %0.2f%%\n'%(100 * tot_trans / max_trans))

    tot_trans = total_transmission(name='Crossbill 2021',
                                   objectives=3,
                                   scan_lenses=2,
                                   mirrors=1,
                                   dichroics=1,
                                   filters=1,
                                   NA1=1.00,
                                   n1=1.33,
                                   NA2=1.00,
                                   n2=1.33,
                                   NA3=1.00,
                                   n3=1.33,
                                   theta_t_deg=45,
                                   verbose=verbose)
    print('Normalised transmission = %0.2f%%\n'%(100 * tot_trans / max_trans))

    tot_trans = total_transmission(name='Lattice 2014',
                                   objectives=1,
                                   scan_lenses=0,
                                   mirrors=1,
                                   dichroics=0,
                                   filters=1,
                                   NA1=1.10,
                                   n1=1.33,
                                   NA2=None,
                                   n2=None,
                                   NA3=None,
                                   n3=None,
                                   theta_t_deg=None,
                                   verbose=verbose)
    print('Normalised transmission = %0.2f%%\n'%(100 * tot_trans / max_trans))

    print('***Coverslip microscopes:***\n')

    tot_trans = total_transmission(name='OPM 2008, built system',
                                   objectives=3,
                                   scan_lenses=0,
                                   mirrors=0,
                                   dichroics=1,
                                   filters=1,
                                   NA1=1.35,
                                   n1=1.51,
                                   NA2=0.85,
                                   n2=1.00,
                                   NA3=0.30,
                                   n3=1.00,
                                   theta_t_deg=32,
                                   verbose=verbose)
    print('Normalised transmission = %0.2f%%\n'%(100 * tot_trans / max_trans))

    tot_trans = total_transmission(name='OPM 2011',
                                   objectives=3,
                                   scan_lenses=0,
                                   mirrors=0,
                                   dichroics=0,
                                   filters=1,
                                   NA1=1.2,
                                   n1=1.33,
                                   NA2=0.95,
                                   n2=1.00,
                                   NA3=0.60,
                                   n3=1.00,
                                   theta_t_deg=32,
                                   verbose=verbose)
    print('Normalised transmission = %0.2f%%\n'%(100 * tot_trans / max_trans))

    tot_trans = total_transmission(name='eSPIM 2019',
                                   objectives=3,
                                   scan_lenses=2,
                                   mirrors=2,
                                   dichroics=1,
                                   filters=1,
                                   NA1=1.27,
                                   n1=1.33,
                                   NA2=0.90,
                                   n2=1.00,
                                   NA3=1.00,
                                   n3=1.33,
                                   theta_t_deg=30,
                                   verbose=verbose)
    print('Normalised transmission = %0.2f%%\n'%(100 * tot_trans / max_trans))

    tot_trans = total_transmission(name='Lattice Zeiss',
                                   objectives=1,
                                   scan_lenses=0,
                                   mirrors=1,
                                   dichroics=0,
                                   filters=1,
                                   NA1=1.00,
                                   n1=1.33,
                                   NA2=None,
                                   n2=None,
                                   NA3=None,
                                   n3=None,
                                   theta_t_deg=None,
                                   verbose=verbose)
    print('Normalised transmission = %0.2f%%\n'%(100 * tot_trans / max_trans))

    tot_trans = total_transmission(name='SOLS 2019',
                                   objectives=3,
                                   scan_lenses=2,
                                   mirrors=2,
                                   dichroics=1,
                                   filters=1,
                                   NA1=1.35,
                                   n1=1.41,
                                   NA2=0.95,
                                   n2=1.00,
                                   NA3=1.00,
                                   n3=1.50,
                                   theta_t_deg=30,
                                   verbose=verbose)
    print('Normalised transmission = %0.2f%%\n'%(100 * tot_trans / max_trans))

    tot_trans = total_transmission(name='SOLS 2 mirror scan, optimum tilt',
                                   objectives=3,
                                   scan_lenses=0,
                                   mirrors=3,
                                   dichroics=1,
                                   filters=1,
                                   NA1=1.35,
                                   n1=1.41,
                                   NA2=0.95,
                                   n2=1.00,
                                   NA3=1.00,
                                   n3=1.50,
                                   theta_t_deg=25,
                                   verbose=verbose)
    print('Normalised transmission = %0.2f%%\n'%(100 * tot_trans / max_trans))

