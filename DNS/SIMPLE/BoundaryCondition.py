 
def VerticalNoSlip(cof_Mat, rows, cols, col_num):

    eq = (rows - 1) * (cols - 1) -1

    for num in range(rows):
        node = cols * num + (col_num - 1)
        cof_Mat[node][node] = 1
        cof_Mat[eq][node] = 0
        
    return cof_Mat
 
def HorizontalNoSlip(cof_Mat, rows, cols, row_num):

    eq = (rows - 1) * (cols - 1) -1

    for num in range(cols):
        node = cols * (row_num - 1) + num
        cof_Mat[node][node] = 1
        cof_Mat[eq][node] = 0
        
    return cof_Mat

def VerticalUltimateVelocity(cof_Mat, rows, cols, col_num):

    u_ult = 5
    eq = (rows - 1) * (cols - 1) -1

    for num in range(rows):
        node = cols * num + (col_num - 1)
        cof_Mat[node][node] = 1
        cof_Mat[eq][node] = u_ult
        
    return cof_Mat

def HorizontalUltimateVelocity(cof_Mat, rows, cols, row_num):

    u_ult = 5
    eq = (rows - 1) * (cols - 1) -1

    for num in range(cols):
        node = cols * (row_num - 1) + num
        cof_Mat[node][node] = 1
        cof_Mat[eq][node] = u_ult
        
    return cof_Mat