from SIMPLE import SIMPLE 
import csv
import numpy as np


p_ps = np.genfromtxt('Results/PreviousStepPressure.csv', delimiter=',')
u_ps = np.genfromtxt('Results/PreviousStepU.csv', delimiter=',')
v_ps = np.genfromtxt('Results/PreviousStepV.csv', delimiter=',')


p_ps = np.array(p_ps).T.tolist()
u_ps = np.array(u_ps).T.tolist()
v_ps = np.array(v_ps).T.tolist()



del_h = 0.0001
del_v = 0.0001
del_t = 0.0001  

time_variation = 1
time_steps = int(time_variation / del_t) + 1 

geo_mat_shape = np.shape(p_ps)
cols = geo_mat_shape[0]
rows = geo_mat_shape[1]

simulation_list = np.empty([rows, cols, 3, time_steps])

flag    = 1
counter = 0

while (flag):

    Ans_feild = SIMPLE(u_ps,v_ps,p_ps,del_h,del_v,del_t)

    for i in range(cols):    
        for j in range(rows):           

            Ans_feild[i][j][0] = u_ps[i][j]
            Ans_feild[i][j][1] = v_ps[i][j]
            Ans_feild[i][j][2] = p_ps[i][j]

            simulation_list[i][j][0][counter] = Ans_feild[i][j][0]
            simulation_list[i][j][1][counter] = Ans_feild[i][j][1]
            simulation_list[i][j][2][counter] = Ans_feild[i][j][2]

    counter += 1

    if counter == time_steps:
        flag = 0

flow_simulation = np.array(simulation_list)
np.savetxt('Results/low_simulation.txt', flow_simulation, fmt='%f')
