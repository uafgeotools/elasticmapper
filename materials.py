import numpy as np
import sys

def reorder(input_vec, output_index):
    output_vec = np.zeros(21)
    for i, j in enumerate(output_index):
        output_vec[j - 1] = input_vec[i]
    return output_vec

def read_material(filename):
    def read_table(filename):
        table = []
        with open(filename, 'r') as file:
            for line in file:
                row = [item for item in line.split()]
                table.append(row)
        return table

    def reorganize(table):
        if "Vestrum" in filename:
            cvec = np.array(   table[0][0:6]
                             + table[1][1:6]
                             + table[2][2:6]
                             + table[3][3:6]
                             + table[4][4:6]
                             + table[5][5:6])

            stdv = np.array(   table[0][6: ]
                             + table[1][7: ]
                             + table[2][8: ]
                             + table[3][9: ]
                             + table[4][10:]
                             + table[5][11:])

            cvec = [float(item) for item in cvec]
            stdv = [float(item) for item in stdv]
            materials = np.array(list(zip(cvec, stdv)))
            order = [i for i in range(1,22)]

        elif any(item in filename for item in ["Aminzadeh", "Brown", "Lokajicek",
                                                    "Militzer"]):
            if "BrownAbramson" in filename:
                materials = []
                for i in range(len(table)):
                    row = []
                    for j in range(1,len(table[0])-1,2):
                        row.append([float(table[i][j]), float(table[i][j+1])])
                    materials.append(row)
            elif "Brown" in filename:
                materials = []
                for i in range(21):
                    row = []
                    for j in range(1,16,2):
                        row.append([float(table[i][j]), float(table[i][j+1])])
                    materials.append(row)
            else:
                materials = []
                for row in table:
                    materials.append([float(item) for item in row[1:]])

            cmn = ["C11", "C12", "C13", "C14", "C15", "C16",
                          "C22", "C23", "C24", "C25", "C26",
                                 "C33", "C34", "C35", "C36",
                                        "C44", "C45", "C46",
                                               "C55", "C56",
                                                      "C66"]

            cij = [row[0] for row in table]
            cij = cij + [item for item in cmn if item not in cij]

            materials = np.array(materials)
            order = [np.argmax(np.argsort(cij) == i)+1 for i in range(21)]

        else:
            sys.exit("Reading logic for the given file is not defined")

        return materials, order

    def main(filename):
        table = read_table(filename)
        return reorganize(table)

    return main(filename)

def get_materials_Tvec(material):

    '''
    1 - TapeTape2021 rotated matrices
    2 - TapeTape2021 principal matrices
    3 - TapeTape2022 elasticity tensors from ES_RefMatricesForCollage.nb
    4 - Test matrices from TapeTape2022
    5 - Closest elasticity tensors to a tensor with trivial symmetry far from the monoclinic class (Values taken from ES_FarFromMONO_n2000.pdf)
    6 - Test matrices from exam
    7 - Mar17 matrix
    '''

    s1 = 1
    s2 = np.sqrt(2)
    s3 = np.sqrt(3)
    s6 = np.sqrt(6)
    s8 = np.sqrt(8)

    if material == "TapeTape2021_rotated":

        TRIV = (1/5) * np.array([ 5, 0, 0, -2,  0, 0,
                                     6, 0,  0, -2, 0,
                                        5,  0,  0, 0,
                                            5,  0, 0,
                                                6, 0,
                                                   6])

        MONO = (1/80) * np.array([ 222*s1, -12*s6, -82*s1,  21*s6, -39*s2,   0*s1,
                                           196*s1,  12*s6,  30*s1,  10*s3, -16*s1,
                                                   222*s1,   9*s6, -51*s2,   0*s1,
                                                           242*s1,  -6*s3, -24*s1,
                                                                   254*s1,  -8*s3,
                                                                           304*s1])

        ORTH = (1/10) * np.array([ 38*s1, 18*s1,  3*s2, -5*s2,   -s6,   -s8,
                                          38*s1,  3*s2,  5*s2,   -s6,   -s8,
                                                 41*s1,  0*s1, -7*s3,  6*s1,
                                                        20*s1,  0*s1,  0*s1,
                                                               27*s1, -2*s3,
                                                                      56*s1])

        TET  = (1/64) * np.array([ 168*s1,   4*s6, -40*s1,   6*s6,   6*s2,   0*s1,
                                           324*s1,  -4*s6, -42*s1, -14*s3,  16*s3,
                                                   168*s1,  -6*s6,  -6*s2,   0*s1,
                                                           233*s1,  35*s3,  -8*s3,
                                                                   163*s1,  -8*s1,
                                                                           352*s1])

        CUBE = (1/36) * np.array([ 52*s1,  4*s1, 16*s1, -6*s1, -2*s3,   0*s1,
                                          64*s1,  4*s1, 12*s1,  4*s3,   0*s1,
                                                 52*s1, -6*s1, -2*s3,   0*s1,
                                                        45*s1,  3*s3,   0*s1,
                                                               39*s1,   0*s1,
                                                                      108*s1])

        TRIG = (1/16) * np.array([ 32-4*s3,   -8*s1,  3*s2, -6*s2, -3*s6,  0*s1,
                                            32+4*s3,  3*s2,  6*s2, -3*s6,  0*s1,
                                                     76*s1, -2*s3, 12*s3,  4*s3,
                                                            24*s1,  6*s1,  0*s1,
                                                                   52*s1,  4*s1,
                                                                          88*s1] )

        XISO = (1/128) * np.array([ 532*s1,  92*s3, 122*s1,  -2*s3,  -60*s1,  48*s1,
                                            348*s1,  -2*s3, 126*s1,  -20*s3,  16*s3,
                                                    349*s1,  31*s3, -126*s1,  24*s1,
                                                            287*s1,  -42*s3,   8*s3,
                                                                     212*s1, -16*s1,
                                                                             704*s1])

        ISO  = np.array([ 1, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0,
                                1, 0, 0, 0,
                                   1, 0, 0,
                                      1, 0,
                                         3])

        return TRIV, MONO, ORTH, TET, CUBE, TRIG, XISO, ISO

    if material == "TapeTape2021_principal":

        TRIV = (1/5) * np.array([ 5, 0, 0, -2,  0, 0,
                                     6, 0,  0, -2, 0,
                                        5,  0,  0, 0,
                                            5,  0, 0,
                                                6, 0,
                                                   6])

        MONO = (1/20) * np.array([ 50+15*s3,      -15,  0,  0,  0,  0,
                                             50-15*s3,  0,  0,  0,  0,
                                                       64,  0,  0, -8,
                                                           76, 12,  0,
                                                               44,  0,
                                                                   76])

        ORTH = (1/5) * np.array([ 10,  0, 0,  0,     0,     0,
                                      15, 0,  0,     0,     0,
                                          5,  0,     0,     0,
                                             28,  2*s3,     2,
                                                    24, -2*s3,
                                                           28])

        TET  = (1/2) * np.array([ 4, 0, 0, 0,  0,  0,
                                     4, 0, 0,  0,  0,
                                        8, 0,  0,  0,
                                           6,  0,  0,
                                              11, -1,
                                                  11])

        CUBE = np.array([ 1, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0,
                                1, 0, 0, 0,
                                   2, 0, 0,
                                      2, 0,
                                         3])

        TRIG = (1/2) * np.array([ 3, 0, s3,  0,  0,  0,
                                     3,  0, s3,  0,  0,
                                         5,  0,  0,  0,
                                             5,  0,  0,
                                                11, -1,
                                                    11])

        XISO = (1/2) * np.array([ 2, 0, 0, 0,  0,  0,
                                     2, 0, 0,  0,  0,
                                        6, 0,  0,  0,
                                           6,  0,  0,
                                              11, -1,
                                                  11])

        ISO  = np.array([ 1, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0,
                                1, 0, 0, 0,
                                   1, 0, 0,
                                      1, 0,
                                         3])

        return TRIV, MONO, ORTH, TET, CUBE, TRIG, XISO, ISO

    elif material == "TapeTape2022_collage":

        TRIV = np.array([ 4.086056425972008000,  0.010701286581823899, -0.476223441320207000, 0.917018260845806500, -0.818564693235128500, -0.340248312771245900,
                                                 4.672938167924996000, -0.376834406989721540, 0.293522876341346550, -0.962574548964882000, -1.612273033526386400,
                                                                        2.510880008024643400, 0.026713327098956616,  0.512096098591872400,  0.694367252861254700,
                                                                                              4.216229915649152000,  0.384331609719862900,  0.173061012583925920,
                                                                                                                     2.986133113121339600,  0.097534355871519480,
                                                                                                                                            2.527762369307859000])

        MONO = np.array([ 1.250000000000000000, -0.433012701892219300,  0                   ,  0                   ,  0                   , 0                   ,
                                                 1.750000000000000200,  0                   ,  0                   ,  0                   , 0                   ,
                                                                        4.635052873539182000, -0.064686449064606410, -0.750293377816586700, 0.155272382921874950,
                                                                                               4.265283351323687000,  0.776026174574774100, 1.081879633193907800,
                                                                                                                      4.243583384280638500, 0.056114869469269024,
                                                                                                                                            4.856080390856493000])

        ORTH = np.array([ 1, 0, 0, 0                   ,  0                   ,  0                   ,
                             3, 0, 0                   ,  0                   ,  0                   ,
                                5, 0                   ,  0                   ,  0                   ,
                                   2.427256325043446200, -0.323382239501930900,  0.207742446456421100,
                                                          1.614722521590946900, -0.140890855884923730,
                                                                                 3.958021153365611800])

        TET  = np.array([ 1, 0, 0, 0, 0,  0,
                             1, 0, 0, 0,  0,
                                3, 0, 0,  0,
                                   6, 0,  0,
                                      5, -1,
                                          5])

        CUBE = np.array([ 1, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0,
                                1, 0, 0, 0,
                                   2, 0, 0,
                                      2, 0,
                                         3])

        TRIG = np.array([ 1.500000000000000000,  0                   , -0.866025400000000000,  0                   ,  0                   ,  0                   ,
                                                 1.500000000000000000,  0                   , -0.866025400000000000,  0                   ,  0                   ,
                                                                        2.500000000000000000,  0                   ,  0                   ,  0                   ,
                                                                                               2.500000000000000000,  0                   ,  0                   ,
                                                                                                                      5.500000000000000000, -0.500000000000000000,
                                                                                                                                             5.500000000000000000])

        XISO = np.array([ 1, 0, 0, 0, 0, 0  ,
                             1, 0, 0, 0, 0  ,
                                2, 0, 0, 0  ,
                                   2, 0, 0  ,
                                      3, 2.5,
                                         4  ])

        ISO  = np.array([ 1, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0,
                                1, 0, 0, 0,
                                   1, 0, 0,
                                      1, 0,
                                         3])

        return TRIV, MONO, ORTH, TET, CUBE, TRIG, XISO, ISO

    elif material == "TapeTape2022_test":

        T1vec = (1/128) * np.array([ 532*s1,  92*s3, 122*s1,  -2*s3,  -60*s1,  48*s1,
                                             348*s1,  -2*s3, 126*s1,  -20*s3,  16*s3,
                                                     349*s1,  31*s3, -126*s1,  24*s1,
                                                             287*s1,  -42*s3,   8*s3,
                                                                      212*s1, -16*s1,
                                                                              704*s1])

        T2vec = (1/64) * np.array([ 168*s1,   4*s6, -40*s1,   6*s6,   6*s2,   0*s1,
                                            324*s1,  -4*s6, -42*s1, -14*s3,  16*s3,
                                                    168*s1,  -6*s6,  -6*s2,   0*s1,
                                                            233*s1,  35*s3,  -8*s3,
                                                                    163*s1,  -8*s1,
                                                                            352*s1])

        T3vec = (1/320) * np.array([ 808*s1,    4*s6, -168*s1,    6*s6,   6*s2,    0*s1,
                                             1572*s1,   -4*s6, -282*s1, -94*s3,   80*s3,
                                                       808*s1,   -6*s6,  -6*s2,    0*s1,
                                                               1057*s1, 139*s3,  -40*s3,
                                                                        779*s1,  -40*s1,
                                                                                1760*s1])

        return T1vec, T2vec, T3vec

    elif material == "FarFromMONO":

        TRIV  = np.array([ 0.31708200, -0.00658517,  0.22786200,  0.21917600,  0.11476200, -0.24482300,
                                        0.07346000, -0.11485800, -0.02070150,  0.02757740, -0.14234200,
                                                     0.42042400,  0.21279300,  0.05286390,  0.12043800,
                                                                  0.24103100,  0.10226500, -0.16281800,
                                                                               0.09627170, -0.14736100,
                                                                                            0.64901000])

        MONO  = np.array([ 0.30453100, -0.04445250,  0.23883800,  0.18417100,  0.09108580, -0.18169800,
                                        0.10199400, -0.16294800, -0.04064210,  0.03111810, -0.14669500,
                                                     0.44764100,  0.19559300,  0.00627925,  0.17902700,
                                                                  0.19616900,  0.07340690, -0.09519940,
                                                                               0.09793330, -0.18417100,
                                                                                            0.64901000])

        ORTH  = np.array([ 0.32819100, -0.01518270,  0.23797500,  0.25159800,  0.03564770, -0.15428700,
                                        0.13682700, -0.06213200, -0.02213080,  0.04212710, -0.21736400,
                                                     0.31050700,  0.22668300,  0.02508750, -0.01148220,
                                                                  0.28432000,  0.04669110, -0.13660700,
                                                                               0.08842370, -0.14548700,
                                                                                            0.64901000])

        TET   = np.array([ 0.29491600, -0.04787780,  0.24597800,  0.21226700,  0.01517680, -0.08824300,
                                        0.16588200, -0.10772200, -0.04524340,  0.05081530, -0.24126800,
                                                     0.34913300,  0.22342700, -0.01561020,  0.06071490,
                                                                  0.24903900,  0.01171260, -0.07189800,
                                                                               0.08929750, -0.14812100,
                                                                                            0.64901000])

        CUBE  = np.array([ 0.21051900,  0.04499160,  0.12916000,  0.14135900,  0.05925590,  0.00000000,
                                        0.32138300, -0.10714800,  0.03043980,  0.17098800,  0.00000000,
                                                     0.25167700,  0.12622700, -0.03896390,  0.00000000,
                                                                  0.18948100,  0.04800570,  0.00000000,
                                                                               0.17520800,  0.00000000,
                                                                                            0.64901000])

        TRIG  = np.array([ 0.29898400,  0.07446800,  0.14747000,  0.25023500,  0.12149500, -0.16298400,
                                        0.13417600, -0.00379364,  0.09698110, -0.00111659, -0.04727210,
                                                     0.22286100,  0.15267000,  0.03240060, -0.10343300,
                                                                  0.30724500,  0.08161600, -0.16330700,
                                                                               0.18500300, -0.06859960,
                                                                                            0.64901000])

        XISO  = np.array([ 0.32052500,  0.05793800,  0.15921400,  0.21961300,  0.12242200, -0.15768600,
                                        0.15036500,  0.02560920,  0.08364730,  0.03772530, -0.04859230,
                                                     0.19542400,  0.16733600,  0.05284360, -0.10746900,
                                                                  0.32719800,  0.07759900, -0.15781400,
                                                                               0.15475600, -0.06907050,
                                                                                            0.64901000])

        ISO   = np.array([ 0.22965400,  0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                        0.22965400,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                                     0.22965400,  0.00000000,  0.00000000,  0.00000000,
                                                                  0.22965400,  0.00000000,  0.00000000,
                                                                               0.22965400,  0.00000000,
                                                                                            0.64901000])

        return TRIV, MONO, ORTH, TET, CUBE, TRIG, XISO, ISO

    elif material == "Test":

        T1vec = (1/16) * np.array([ 32*s1,  0*s1, -16*s1,  0*s1,  0*s1,   0*s1,
                                           60*s1,   0*s1, -6*s1, -2*s3, -40*s3,
                                                   32*s1,  0*s1,  0*s1,   0*s1,
                                                          43*s1,  9*s3,  20*s3,
                                                                 25*s1,  20*s1,
                                                                         32*s1])

        T2vec = np.array([ 1, 0, 0, 0, 0, 0,
                              1, 0, 0, 0, 0,
                                 3, 0, 0, 0,
                                    3, 0, 0,
                                       4, 5,
                                          2])

        T3vec = np.array([ 1, 0, 0, 0, 0, 0,
                              1, 0, 0, 0, 0,
                                 5, 0, 0, 0,
                                    3, 0, 0,
                                       4, 5,
                                          2])

        return T1vec, T2vec, T3vec

    elif material == "Mar17":

        Tvec = np.array([ 206.2,  -67.4,   17.1,  -15.9,  166.6,  -19.5,
                                  352.9,   32.2,    3.3,   61.6,  -41.6,
                                          343.8,    0.4,   73.6,   -4.2,
                                                  270.0,   88.8,   13.3,
                                                          287.1,   30.7,
                                                                   74.0])

        return Tvec

    else:

        sys.exit("Error: requested set index invalid")

def get_materials_Cvec(material):

    '''
    1 - Brown2016 matrices
    2 - Igel1995 matrices
    '''

    if material == "BrownAbramson":

        filename = '../data/BrownAbramson2016_Table3.txt'
        materials, order = read_material(filename)

        output_index = order[:len(materials)]

        cvec_1 = reorder(materials[:, 0, 0], output_index)
        cvec_2 = reorder(materials[:, 1, 0], output_index)
        cvec_3 = reorder(materials[:, 2, 0], output_index)
        cvec_4 = reorder(materials[:, 3, 0], output_index)
        cvec_5 = reorder(materials[:, 4, 0], output_index)
        cvec_6 = reorder(materials[:, 5, 0], output_index)
        cvec_7 = reorder(materials[:, 6, 0], output_index)
        cvec_8 = reorder(materials[:, 7, 0], output_index)
        cvec_9 = reorder(materials[:, 8, 0], output_index)

        stdv_1 = 1/2 * reorder(materials[:, 0, 1], output_index)
        stdv_2 = 1/2 * reorder(materials[:, 1, 1], output_index)
        stdv_3 = 1/2 * reorder(materials[:, 2, 1], output_index)
        stdv_4 = 1/2 * reorder(materials[:, 3, 1], output_index)
        stdv_5 = 1/2 * reorder(materials[:, 4, 1], output_index)
        stdv_6 = 1/2 * reorder(materials[:, 5, 1], output_index)
        stdv_7 = 1/2 * reorder(materials[:, 6, 1], output_index)
        stdv_8 = 1/2 * reorder(materials[:, 7, 1], output_index)
        stdv_9 = 1/2 * reorder(materials[:, 8, 1], output_index)

        return (cvec_1, stdv_1, cvec_2, stdv_2, cvec_3, stdv_3, cvec_4, stdv_4,
                cvec_5, stdv_5, cvec_6, stdv_6, cvec_7, stdv_7, cvec_8, stdv_8,
                cvec_9, stdv_9)

    elif material == "Brown":

        filename = '../data/Brown2016_Table2.txt'
        An, output_index = read_material(filename)

        An0  = reorder(An[:, 0, 0], output_index)
        An25 = reorder(An[:, 1, 0], output_index)
        An37 = reorder(An[:, 2, 0], output_index)
        An48 = reorder(An[:, 3, 0], output_index)
        An60 = reorder(An[:, 4, 0], output_index)
        An67 = reorder(An[:, 5, 0], output_index)
        An78 = reorder(An[:, 6, 0], output_index)
        An96 = reorder(An[:, 7, 0], output_index)

        An0_stdv  = 1/2 * reorder(An[:, 0, 1], output_index)
        An25_stdv = 1/2 * reorder(An[:, 1, 1], output_index)
        An37_stdv = 1/2 * reorder(An[:, 2, 1], output_index)
        An48_stdv = 1/2 * reorder(An[:, 3, 1], output_index)
        An60_stdv = 1/2 * reorder(An[:, 4, 1], output_index)
        An67_stdv = 1/2 * reorder(An[:, 5, 1], output_index)
        An78_stdv = 1/2 * reorder(An[:, 6, 1], output_index)
        An96_stdv = 1/2 * reorder(An[:, 7, 1], output_index)

        return (An0 , An0_stdv , An25, An25_stdv, An37, An37_stdv, An48,
                An48_stdv, An60, An60_stdv, An67, An67_stdv, An78, An78_stdv,
                An96, An96_stdv)

    elif material == "Igel":

        cvec = np.array([10.0,  3.50, 2.50, -5.000,  0.10,  0.300,
                                8.00, 1.50,  0.200, -0.10, -0.150,
                                      6.00,  1.000,  0.40,  0.240,
                                             5.000,  0.35,  0.525,
                                                     4.00, -1.000,
                                                            3.000])

        return cvec

    elif material == "Aminzadeh":

        filename_1 = '../data/Aminzadeh2022_Table3.txt'
        BUK, output_index_1 = read_material(filename_1)

        BUK_01  = reorder(BUK[:,0], output_index_1)
        BUK_2   = reorder(BUK[:,1], output_index_1)
        BUK_5   = reorder(BUK[:,2], output_index_1)
        BUK_10  = reorder(BUK[:,3], output_index_1)
        BUK_20  = reorder(BUK[:,4], output_index_1)
        BUK_50  = reorder(BUK[:,5], output_index_1)
        BUK_80  = reorder(BUK[:,6], output_index_1)
        BUK_100 = reorder(BUK[:,7], output_index_1)

        filename_2 = '../data/Aminzadeh2022_Table4.txt'
        GRM, output_index_2 = read_material(filename_2)

        GRM_01  = reorder(GRM[:, 0], output_index_2)
        GRM_5   = reorder(GRM[:, 1], output_index_2)
        GRM_10  = reorder(GRM[:, 2], output_index_2)
        GRM_15  = reorder(GRM[:, 3], output_index_2)
        GRM_20  = reorder(GRM[:, 4], output_index_2)
        GRM_50  = reorder(GRM[:, 5], output_index_2)
        GRM_80  = reorder(GRM[:, 6], output_index_2)
        GRM_100 = reorder(GRM[:, 7], output_index_2)

        return (BUK_01, BUK_2, BUK_5, BUK_10, BUK_20, BUK_50, BUK_80, BUK_100,
                GRM_01, GRM_5, GRM_10, GRM_15, GRM_20, GRM_50, GRM_80, GRM_100)

    elif material == "Lokajicek":

        filename = '../data/Lokajicek2021_supp_5MPa.txt'
        materials, output_index_1 = read_material(filename)

        WG2D = materials[:, [3,6,8,9]]

        WG2D_100 = reorder(WG2D[:, 0], output_index_1)
        WG2D_200 = reorder(WG2D[:, 1], output_index_1)
        WG2D_400 = reorder(WG2D[:, 2], output_index_1)
        WG2D_600 = reorder(WG2D[:, 3], output_index_1)

        filename = '../data/Lokajicek2021_Table3.txt'
        materials, output_index_2 = read_material(filename)

        WG100 = materials[:, [4, 5]]
        WG600 = materials[:, [6, 7]]

        WG100_01  = reorder(WG100[:, 0], output_index_2)
        WG100_400 = reorder(WG100[:, 1], output_index_2)
        WG600_01  = reorder(WG600[:, 0], output_index_2)
        WG600_400 = reorder(WG600[:, 1], output_index_2)

        return (WG2D_100, WG2D_200, WG2D_400, WG2D_600, WG100_01, WG100_400,
                WG600_01, WG600_400)

    elif material == "Vestrum":

        filename = "../data/Vestrum1996_Table2.txt"
        material, output_index = read_material(filename)
        cvec = reorder(material[:, 0], output_index)
        stdv = reorder(material[:, 1], output_index)

        return cvec, stdv

    else:

        sys.exit("Error: requested set index invalid")
