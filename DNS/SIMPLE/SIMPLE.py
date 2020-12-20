import numpy as np


def SIMPLE(u_ps,v_ps,p_ps,del_h,del_v,del_t):

    from GaussianElimination import gauss
    import BoundaryCondition

    geo_mat_shape = np.shape(p_ps)
    cols = geo_mat_shape[0]
    rows = geo_mat_shape[1]

    abs_Conv   = 0

    #Define an Answer Array in 3D to be the Return Parameter
    Ans = [[[0 for j in range(rows)] for i in range(cols)]for k in range(3)]

    while(abs_Conv == False):

        Au   = [[0 for j in range(rows)] for i in range(cols)]
        Au_n = [[0 for j in range(rows)] for i in range(cols)]
        Au_s = [[0 for j in range(rows)] for i in range(cols)]
        Au_w = [[0 for j in range(rows)] for i in range(cols)]
        Au_e = [[0 for j in range(rows)] for i in range(cols)]
        b_u  = [[0 for j in range(rows)] for i in range(cols)]

        Av   = [[0 for j in range(rows)] for i in range(cols)]
        Av_n = [[0 for j in range(rows)] for i in range(cols)]
        Av_s = [[0 for j in range(rows)] for i in range(cols)]
        Av_w = [[0 for j in range(rows)] for i in range(cols)]
        Av_e = [[0 for j in range(rows)] for i in range(cols)]
        b_v  = [[0 for j in range(rows)] for i in range(cols)]

        Ap   = [[0 for j in range(rows)] for i in range(cols)]
        Ap_n = [[0 for j in range(rows)] for i in range(cols)]
        Ap_s = [[0 for j in range(rows)] for i in range(cols)]
        Ap_w = [[0 for j in range(rows)] for i in range(cols)]
        Ap_e = [[0 for j in range(rows)] for i in range(cols)]
        b_p  = [[0 for j in range(rows)] for i in range(cols)]


        node_num = rows * cols #!!!Just for Rectangular Geometry!!!

        cof_Mat_rows = (node_num)
        cof_Mat_cols = (node_num + 1)

        u_cof_Mat = [[0 for i in range(cof_Mat_cols)] for j in range(cof_Mat_rows)]
        v_cof_Mat = [[0 for i in range(cof_Mat_cols)] for j in range(cof_Mat_rows)]
        p_cof_Mat = [[0 for i in range(cof_Mat_cols)] for j in range(cof_Mat_rows)]# delta_P!!!

        
        u_s     = [[0 for j in range(rows)] for i in range(cols)]
        v_s     = [[0 for j in range(rows)] for i in range(cols)]
        u_c     = [[0 for j in range(rows)] for i in range(cols)]
        v_c     = [[0 for j in range(rows)] for i in range(cols)]
        delta_P = [[0 for j in range(rows)] for i in range(cols)]



        conv_elmt1 = 1  #Define Elements for Convergence
        conv_elmt2 = 0        
        absErr = 0.00001 #Desired Error for dela_P
        varErr = 0.00001 #Desired Error for Variation
        varVal = 0       #To calculate variation Error Value



        for i in range(0, (cols - 1)):
            for j in range(0, (rows - 1)):

                G_node_num = cols * j + i


                #del_t = 0.00 #Time Step
                #del_h = 0.00 #Horizontal Step
                #del_v = 0.00 #Vertical   Step

                #Calculation Some Constants due to a More Efficient use of Memory
                #velocity = ((u_ps[i][j]) * (u_ps[i][j]) + (v_ps[i][j]) *(v_ps[i][j])) ** 0.5
                #rho = 0.01
                #L = 0.01
                #Mu = 0.01
                Re = 100 #(rho * velocity * L) / Mu

                h_vr = del_h / (del_v * Re)
                v_hr = del_v / (del_h * Re)

                hv_t = (del_h * del_v) / del_t

                #Horizental Momentum coefficients

                #u_ps represents horizontal direction velocity of previous step
                #v_ps represents vertical   direction velocity of previous step

                Au[i][j] = ( 0.25 *( del_v * ( (u_ps[i][j] + u_ps[i+1][j]) - (u_ps[i][j] + u_ps[i-1][j]) ) + \
                    del_h * ( (v_ps[i][j] + v_ps[i+1][j]) - (v_ps[i][j] +v_ps[i-1][j]) )) \
                     + 2 * (v_hr * h_vr) + hv_t)  

                Au_n[i][j] =   0.25 * del_h * (v_ps[i][j] + v_ps[i+1][j]) - v_hr
                Au_s[i][j] =  -0.25 * del_h * (v_ps[i][j-1] + v_ps[i+1][j-1]) - v_hr
                Au_w[i][j] =  -0.25 * del_v * (u_ps[i][j] + u_ps[i-1][j]) -h_vr
                Au_e[i][j] =   0.25 * del_v * (u_ps[i][j] + u_ps[i+1][j]) - h_vr

                b_u[i][j]  = hv_t * u_ps[i][j]



                u_cof_Mat[G_node_num][G_node_num] = Au[i][j]

                u_cof_Mat[G_node_num - cols][G_node_num] = Au_n[i][j]
                u_cof_Mat[G_node_num + cols][G_node_num] = Au_s[i][j]
                u_cof_Mat[G_node_num - 1][G_node_num]    = Au_w[i][j]
                u_cof_Mat[G_node_num + 1][G_node_num]    = Au_e[i][j]

                u_cof_Mat[G_node_num][node_num] = -b_u[i][j] - del_h * (p_ps[i+1][j] - p_ps[i][j])

                #Vertical   Momentum coefficients

                Av[i][j] = ( 0.25 *( del_v * ( (u_ps[i-1][j] + u_ps[i-1][j+1]) - (u_ps[i][j] + u_ps[i][j+1]) ) + \
                    del_h * ( (v_ps[i][j] + v_ps[i][j+1]) - (v_ps[i][j-1] + v_ps[i][j]) )) + \
                        2 * (v_hr * h_vr) + hv_t)

                Av_n[i][j] =   0.25 * del_h * (v_ps[i][j]   + v_ps[i][j+1])   - v_hr
                Av_s[i][j] =  -0.25 * del_h * (v_ps[i][j-1] + v_ps[i][j])     - v_hr  
                Av_w[i][j] =  -0.25 * del_v * (u_ps[i-1][j] + u_ps[i-1][j+1]) -h_vr
                Av_e[i][j] =   0.25 * del_v * (u_ps[i][j]   + u_ps[i][j+1])   - h_vr

                b_v[i][j]  = hv_t * v_ps[i][j]



                v_cof_Mat[G_node_num][G_node_num] = Av_n[i][j]

                v_cof_Mat[G_node_num - cols][G_node_num] = Av_n[i][j]
                v_cof_Mat[G_node_num + cols][G_node_num] = Av_s[i][j]
                v_cof_Mat[G_node_num - 1][G_node_num]    = Av_w[i][j]
                v_cof_Mat[G_node_num + 1][G_node_num]    = Av_e[i][j]

                v_cof_Mat[G_node_num][(node_num)] = -b_v[i][j] - del_v * (p_ps[i][j+1] - p_ps[i][j])


        u_cof_Mat = BoundaryCondition.HorizontalNoSlip(u_cof_Mat, rows, cols, 0)
        u_cof_Mat = BoundaryCondition.HorizontalUltimateVelocity(u_cof_Mat, rows, cols, (rows - 1))        
        v_cof_Mat = BoundaryCondition.VerticalNoSlip(v_cof_Mat, rows, cols, 0)
        v_cof_Mat = BoundaryCondition.VerticalNoSlip(v_cof_Mat, rows, cols, (cols - 1))
        
        u_s_G = gauss(u_cof_Mat)   #Solves Matrix for a 1D Global Velocity Array
        v_s_G = gauss(v_cof_Mat)


        for i in range(cols):    #To have a 2D Array(Matrix Form) of Velocities
            for j in range(rows):

                u_s[i][j] = u_s_G[(i*node_num + j)]
                v_s[i][j] = v_s_G[(i*node_num + j)]
    
      


        for i in range(cols):
            for j in range(rows):


                #Pressure Correction coefficients

                Ap[i][j] = del_t /((1 + (del_t * Av[i][j]) / (del_h * del_v) ) * del_v) +\
                   del_t /((1 + (del_t * Av[i][j-1]) / (del_h * del_v) ) * del_v) +\
                       del_t /((1 + (del_t * Au[i-1][j]) / (del_h * del_v) ) * del_h) +\
                          del_t /((1 + (del_t * Au[i][j]) / (del_h * del_v) ) * del_h)

                Ap_n[i][j] = -del_t /((1 + (del_t * Av[i][j])   / (del_h * del_v) ) * del_v)
                Ap_s[i][j] = -del_t /((1 + (del_t * Av[i][j-1]) / (del_h * del_v) ) * del_v)
                Ap_w[i][j] = -del_t /((1 + (del_t * Au[i-1][j]) / (del_h * del_v) ) * del_h)
                Ap_e[i][j] = -del_t /((1 + (del_t * Au[i][j])   / (del_h * del_v) ) * del_h)

                b_p[i][j] = -(u_s[i][j] - u_s[i-1][j]) * del_v -\
                    (v_s[i][j] - v_s[i][j-1])



                p_cof_Mat[G_node_num][G_node_num] = Ap_n[i][j]

                p_cof_Mat[G_node_num][G_node_num - cols] = Ap_n[i][j]
                p_cof_Mat[G_node_num][G_node_num + cols] = Ap_s[i][j]
                p_cof_Mat[G_node_num][G_node_num - 1]    = Ap_w[i][j]
                p_cof_Mat[G_node_num][G_node_num + 1]    = Ap_e[i][j]

                p_cof_Mat[G_node_num][(node_num + 1)] = -b_p[i][j]


                if abs(b_p[i][j])<absErr:
                    conv_elmt1 = 0
        
                varVal += ((b_p[i][j]) ** 2)


        delta_P_G = gauss(p_cof_Mat)

        for i in range(cols):    #To have a 2D Array(Matrix Form) of the Delta_P
            for j in range(rows):

                delta_P[i][j] = delta_P_G[(i*node_num + j)]
        
  
        for i in range(cols):    #To Estmaite (u_c, v_c) and calculate New Velocities.
            for j in range(rows):

                u_c[i][j] = (del_t /((1 + (del_t * Au[i][j]) / (del_h * del_v) ) * del_h)) * (delta_P[i][j] - delta_P[i+1][j])
                v_c[i][j] = (del_t /((1 + (del_t * Av[i][j]) / (del_h * del_v) ) * del_v)) * (delta_P[i][j] - delta_P[i][j+1])

                u_ps[i][j] = u_s[i][j] + u_c[i][j]   
                v_ps[i][j] = v_s[i][j] + v_c[i][j]

                p_ps[i][j] += delta_P[i][j]   


        if varVal<varErr:
            conv_elmt2 = 1

        abs_Conv = conv_elmt1 * conv_elmt2 #if Converged abs_Conv would be equal to 1



    for i in range(cols):    
        for j in range(rows):
            Ans[i][j][0] = u_ps[i][j]
            Ans[i][j][1] = v_ps[i][j]
            Ans[i][j][2] = p_ps[i][j]


             



    return Ans

    