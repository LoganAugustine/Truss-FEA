def show_truss(q_sol, elem_num, elems_info, nodes_info):
    ''' function entries:
        1. q_sol: LIST of solutions to global nodal displacements
             - format: [qx_node1, qy_node1, qx_node2, qy_node2, ...,qx_nodeN, qy_nodeN],
               where; - qx_nodeN represents displacement of node N in x-direction
                      - qy_nodeN represents displacement of node N in y-direction

        2. elem_num: number of elements in truss
             - format: An integer number

        3. elems_info: element connectivity information
             - format: [[elem1, elem1_node1, elem1_node2], [elem2, elem2_node1, elem2_node2], ...,[elemN, elemN_node1, elemN_node2]]
                 where; - elemN represents the element number N (an integer)
                        - elemN_node1 represents node1 of element N (an integer)
                        - elemN_node2 represents node2 of element N (an integer)

        4. nodes_info: node global coordinates
             - format: [[node1, x1_coord, y1_coord], [node2, x2_coord, y2_coord], ...,[nodeN, xN_coord, yN_coord]]
                 where; - nodeN represents node number N (an integer)
                        - xN_coord represents x location of node N
                        - yN_coord represents y location of node N

        For Example: JDW Example - Chapter 4, page 14;
              elems_info = [[1, 1, 2], [2, 1, 4], [3, 3, 2], [4, 4, 2], [5, 3, 4]]
              nodes_info = [[1, 0, 10], [2, 10, 10], [3, 0, 0], [4, 10, 0]]
              q_sol = [0, 0, -0.250000000000000, 0.957106781186548, 0, 0, 0, 0.707106781186548]'''

    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    # -------------------------------------------------------------------------
    figt, tplt = plt.subplots(1, 1)
    for e in range(elem_num):
        e_info = elems_info[e]
        e_nodes = [e_info[1], e_info[2]]
        # ---------------------------------------------------------------------
        # getting element dofs
        e_dofs = []
        for i in range(len(e_nodes)):
            e_dofs.append(2 * e_nodes[i] - 1)
            e_dofs.append(2 * e_nodes[i])
        # ---------------------------------------------------------------------
        # computing node length
        node1_info, node2_info = nodes_info[e_nodes[0] - 1], nodes_info[e_nodes[1] - 1]
        node1_x, node1_y, node2_x, node2_y = node1_info[1], node1_info[2], node2_info[1], node2_info[2]
        le_final = ((node2_x - node1_x) ** 2 + (node2_y - node1_y) ** 2) ** 0.5
        le_final = float(le_final)
        # ---------------------------------------------------------------------
        # getting element nodal displacements
        qe_final = []
        for dof in e_dofs:
            qe_final.append(q_sol[dof - 1])
            # print(dof, q_sol[dof - 1])
        # ---------------------------------------------------------------------
        # deformed_coord = initial_coord + deformation, qe
        x_undef_range, y_undef_range = [node1_x, node2_x], [node1_y, node2_y]
        x_def_range, y_def_range = [node1_x + float(qe_final[0]), node2_x + float(qe_final[2])], [
            node1_y + float(qe_final[1]), node2_y + float(qe_final[3])]
        # ---------------------------------------------------------------------
        # Plots
        if e != elem_num - 1:
            tplt.plot(x_undef_range, y_undef_range, 'g--')
            tplt.plot(x_def_range, y_def_range, color='black')
        else:
            tplt.plot(x_undef_range, y_undef_range, 'g--', label='Undeformed')
            tplt.plot(x_def_range, y_def_range, color='black', label='Deformed')

    tplt.set(xlabel='x, mm', ylabel='y, mm')
    tplt.set_title('Truss Deformation')
    tplt.grid()
    tplt.legend()
    plt.show()
    return figt
