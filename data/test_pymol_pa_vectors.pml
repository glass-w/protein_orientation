
                from cgo import *
                axis1=  [ BEGIN, LINES, COLOR, 1.0, 0.0, 0.0,                 VERTEX,    0.000,    0.000,    0.000, VERTEX,   60.000,    0.000,    0.000, END ]
                axis2=  [ BEGIN, LINES, COLOR, 0.0, 1.0, 0.0,                 VERTEX,    0.000,    0.000,    0.000, VERTEX,    0.000,   40.000,    0.000, END ]
                axis3=  [ BEGIN, LINES, COLOR, 0.0, 0.0, 1.0,                 VERTEX,    0.000,    0.000,    0.000, VERTEX,    0.000,    0.000,   20.000, END ]
                cmd.load_cgo(axis1, 'axis1')
                cmd.load_cgo(axis2, 'axis2')
                cmd.load_cgo(axis3, 'axis3')
                cmd.set('cgo_line_width', 4)
                