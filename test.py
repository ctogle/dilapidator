import dilap.construct as dlc

def cube():
    l = 2.4
    dlc.build(dlc.cube(l).translate_x(3).translate_y(0.5).translate_z(5))
    
cube()

