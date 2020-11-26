scale_factor = 20


def vis_axes(vis, axes_data, center, name):

    '''
    Visualise the principal axes being used for the calculation.

    vis : the format to write the axes out as (either 'vmd' or 'pymol')
    axes_data :  the array containing the three principal axes
    center : the geometric center of the users selection
    name : the name of the .pdb file (based on the supplied .xtc file)
    '''

    axis1 = axes_data[0]
    axis2 = axes_data[1]
    axis3 = axes_data[2]

    if vis == 'vmd':

        output = open(str(name) + '_pa_vectors.pdb', 'a')

        for i in range(0, (3 * scale_factor)):
            tmp = "ATOM    {0:3d}  CA  ALA A {1:3d}    {2:8.3f}{3:8.3f}{4:8.3f}\
                  1.00  0.00\n".format(i, i, center[0] + (axis1[0] * i),
                                       center[1] + (axis1[1] * i),
                                       center[2] + (axis1[2] * i))
            output.write(tmp)

        output.write("TER\n")

        for j in range(0, (2 * scale_factor)):
            tmp2 = "ATOM    {0:3d}  CA  ALA B {1:3d}    {2:8.3f}{3:8.3f}{4:8.3f}\
                  1.00  0.00\n".format(j, j, center[0] + (axis2[0] * j),
                                       center[1] + (axis2[1] * j),
                                       center[2] + (axis2[2] * j))
            output.write(tmp2)

        output.write("TER\n")

        for k in range(0, (1 * scale_factor)):
            tmp3 = "ATOM    {0:3d}  CA  ALA C {1:3d}    {2:8.3f}{3:8.3f}{4:8.3f}\
                  1.00  0.00\n".format(k, k, center[0] + (axis3[0] * k),
                                       center[1] + (axis3[1] * k),
                                       center[2] + (axis3[2] * k))
            output.write(tmp3)

        output.write("TER\nENDMDL\n")

        output.close()

    elif vis == 'pymol':

        # --------------------------------------------------------------------------
        # center axes to the geometric center of the molecule
        # and rescale them by order of eigen values
        # --------------------------------------------------------------------------

        # the large vector is the first principal axis
        point1 = 3 * scale_factor * axis1 + center
        # the medium vector is the second principal axis
        point2 = 2 * scale_factor * axis2 + center
        # the small vector is the third principal axis
        point3 = 1 * scale_factor * axis3 + center

        # pymol_name = pdb_name.replace(".pdb", "_axes.pml")

        pymol_name = (name + "_pa_vectors.pml")
        with open(pymol_name, "w") as pymol_file:
            pymol_file.write(
                """
                from cgo import *
                axis1=  [ BEGIN, LINES, COLOR, 1.0, 0.0, 0.0, \
                VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
                axis2=  [ BEGIN, LINES, COLOR, 0.0, 1.0, 0.0, \
                VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
                axis3=  [ BEGIN, LINES, COLOR, 0.0, 0.0, 1.0, \
                VERTEX, %8.3f, %8.3f, %8.3f, VERTEX, %8.3f, %8.3f, %8.3f, END ]
                cmd.load_cgo(axis1, 'axis1')
                cmd.load_cgo(axis2, 'axis2')
                cmd.load_cgo(axis3, 'axis3')
                cmd.set('cgo_line_width', 4)
                """ % (
                    center[0], center[1], center[2], point1[0],
                    point1[1], point1[2],
                    center[0], center[1], center[2], point2[0],
                    point2[1], point2[2],
                    center[0], center[1], center[2], point3[0],
                    point3[1], point3[2]))

        # --------------------------------------------------------------------------
        # create .pml script for nice rendering in Pymol
        # output usage
        # --------------------------------------------------------------------------
        print("\nFirst principal axis (in red)")
        # print("coordinates: ", axis1)
        # print("eigen value: ", eval1)

        print("\nSecond principal axis (in green)")
        # print("coordinates:", axis2)
        # print("eigen value:", eval2)

        print("\nThird principal axis (in blue)")
