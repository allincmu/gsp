import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee
import random
import cProfile
import os
import sys
import datetime
import sympy

folder = ""
timestamp = ""

NODES = 100
PROB = 0.03
SHOW_GRAPHS = True


numpy.set_printoptions(threshold=numpy.inf, suppress=True, precision=2, linewidth=NODES*NODES)

###############################################
##               GSP Interface               ##
###############################################

def redirect_to_file(filename):
    original = sys.stdout
    sys.stdout = open(folder + filename + timestamp + ".txt", 'w')
    return original

def make_dir(direct, time=False, set_folder = False):
    global timestamp
    
    current_date_and_time = datetime.datetime.now()
    date_time = current_date_and_time.strftime("%m%d%y-%H%M")
    date_only = current_date_and_time.strftime("%m%d%y")
    timestamp = str(date_time).replace(" ", "_")

    path = os.getcwd()
    if time:
        name = path + "/" + date_only + "/" + direct + "__" + timestamp
    else:
        name = path + "/" + "/" + direct
    try:
        os.makedirs(name)
    except OSError:
        print ("Creation of the directory %s failed" % name)
    else:
        print ("Successfully created the directory %s" % name)

    if set_folder:
        global folder
        folder = date_only + '/' + direct + "__" +timestamp+"/"
    return direct+timestamp

def gen_conn_graph(nodes = NODES, prob = PROB):
    graph = nx.binomial_graph(NODES, PROB)
    while( not nx.is_connected(graph)):
        graph = nx.binomial_graph(NODES, PROB)
    return graph

def draw_graph(graph, name, show = True):
    print(f"({name}) # edges: ", len(graph.edges()))
    if (show):
        try:
            colors = [(weight) for u, v, weight in graph.edges.data("weight", default=1)]
            # print(colors)
            print('max_edge_weight min_edge_weight', max(colors), min(colors))
            vmin = -1.5
            vmax = 1.5
        except:
            colors = [10] * len(graph.edges())
        cmap = plt.cm.nipy_spectral
        options = {
            "node_color": "#A0CBE2",
            "edge_color": colors,
            "edge_vmin": vmin,
            "edge_vmax": vmax,
            "width": 1,
            "alpha": 0.8,
            "edge_cmap": cmap, #jet
            "with_labels": False,
            "node_size": 20,
        }
        nx.draw_circular(graph, **options)
        sm = plt.cm.ScalarMappable(cmap=cmap, 
                                   norm=plt.Normalize(vmin = vmin, vmax=vmax))
        sm._A = []
        plt.colorbar(sm)

        
def get_evd(adj_mtx):

    # get eigenvalues and vectors and calculate GFT
    eigenvalues, eigenvectors = numpy.linalg.eig(adj_mtx)
    GFT = numpy.linalg.inv(eigenvectors)
    GFT_INV = eigenvectors

    # diagonalize eigenvalues
    h = numpy.diag(eigenvalues)

    # print("h", h)
    # print("GFT", GFT)
    
    return GFT, h, GFT_INV

def get_gft_graph(graph):
    # get adjacency matrix
    adj_mtx = nx.to_numpy_array(graph)
    # print(adj_mtx)

    # get diagonalization of adj_mtx
    GFT, h, GFT_INV = get_evd(adj_mtx)

    # get forier transfromed graph
    ft_mtx = numpy.matmul(numpy.matmul(GFT, numpy.conjugate(h)), GFT_INV)

    # ensure M is symetric
    if(not (ft_mtx.transpose() == ft_mtx).all()):
        ft_mtx = 1/2 * (ft_mtx + ft_mtx.transpose())
        assert((ft_mtx.transpose() == ft_mtx).all())


    # get new graph
    gft_graph = nx.from_numpy_matrix(ft_mtx)

    return gft_graph, GFT, h, GFT_INV

def threshold(graph, threshold):

    # get adjacency matrix
    adj_mtx = nx.to_numpy_array(graph)

    for i in range(len(adj_mtx)):
        for j in range(len(adj_mtx[i])):
            if (abs(adj_mtx[i][j]) < threshold):
                adj_mtx[i][j] = 0

    if(not (adj_mtx.transpose() == adj_mtx).all()):
        adj_mtx = 1/2 * (adj_mtx + adj_mtx.transpose())
        # print(adj_mtx)
        assert((adj_mtx.transpose() == adj_mtx).all())
    

    new_graph = nx.from_numpy_array(adj_mtx)

    if not nx.is_connected(new_graph):
        print("Graph no longer connected")
        print_graphs()
        quit()

    edges_removed = len(graph.edges) - len(new_graph.edges)
    print(f"Edges Removed: {edges_removed}")

    return new_graph

# returns sorted eigenval and eigenvec
def get_evd_of_graph(graph):
    mtx = nx.to_numpy_array(graph)
    if(not (mtx.transpose() == mtx).all()):
        mtx = 1/2 * (mtx + mtx.transpose())
        # print(mtx)
        assert((mtx.transpose() == mtx).all())
    eigenvalues, eigenvectors = numpy.linalg.eig(mtx)

    #### COMMENTED OUT SORTING PART ####

    # if (eigenvalues.max() != 0):
    #     eigenvalues = numpy.abs(1-eigenvalues / numpy.abs(eigenvalues.max()))
    # else:
    #     numpy.abs(eigenvalues)
    
    # eigenvectors = numpy.abs(eigenvectors)

    # eigenvalues, eigenvectors = eigenvalues.tolist(), eigenvectors.tolist()

    # for i in range(10):
    #     a = random.randint(0, NODES-1)
    #     b = random.randint(0, NODES-1)

        
    #     eigen_temp = (eigenvalues[a], eigenvectors[a])
    #     eigenvalues[a], eigenvectors[a] = eigenvalues[b], eigenvectors[b]
    #     eigenvalues[b], eigenvectors[b] = eigen_temp
        
    #     # print(eigenvalues)
    #     # eigenvalues = abs(eigenvalues)
    #     # for i in range(len(eigenvectors)):
    #     #     eigenvectors[i] = abs(eigenvectors[i])


    # # eigenvalues, eigenvectors = zip(*sorted(zip(eigenvalues, eigenvectors)))
    return list(eigenvalues), list(eigenvectors)

def almost_equal(A, B):
    return abs(A-B) < 1e-6

def get_tresholding_error(V_orig, V_th, Thresholding_string):
    plt.close()
    V_orig = numpy.array(V_orig)
    V_orig_tpose = numpy.transpose(V_orig)
    V_th = numpy.array(V_th)
    VthTP = numpy.transpose(V_th)

    # Find Vth' * V
    VthTPxV = numpy.matmul(VthTP, V_orig)

    # Find max val in each row
    max_val = [0] * NODES
    for i in range(len(VthTPxV)):
        max_val[i] = max(abs(VthTPxV[i]))

    # plot max val in each row
    plt.figure()
    plt.scatter(range(NODES), max_val)
    plt.axis([0, NODES, -1.1, 1.1])
    plt.rc('text', usetex=True)
    plt.title("Max value in each row of " + r"$V_{th}^T \cdot V$ " +
              f"(T = {Thresholding_string})", fontsize = 16)
    VthTPxV_copy = VthTPxV.copy()

    # Get error off identity and plot
    err_austin, err_cm = calc_identity_err(VthTPxV_copy)
    plt.text(NODES/2, -0.95, f"Error off Identity (Austin):\n {err_austin:.4f}")
    plt.text(NODES/2, -0.75, f"Error off Identity (Cuthill-Mckee):\n {err_cm:.4f}")
    plt.xlabel("Row", fontsize = 16)
    plt.ylabel("Max Value", fontsize = 16)
    plt.subplots_adjust(left=0.15)
    plt.savefig(f'{folder}Max_Val_Plots/{timestamp}_NA-{NODES}_MaxVal_{Thresholding_string.replace(".", "_")}.png')
    plt.rc('text', usetex=False)



    ###### Find angle Histogram ######
    angles = [0] * NODES * NODES
    VthTPxV = numpy.divide(VthTPxV, numpy.max(numpy.abs(VthTPxV))) #find angles
    try:
        for i in range(len(VthTPxV)):
            for j in range(len(VthTPxV[i])):
                angles[i*len(VthTPxV) + j] = numpy.degrees(numpy.arccos(VthTPxV[i][j]))
    except:
        print(V_th)
        raise(TypeError)
        exit()

    # Plot histogram
    plt.figure()
    plt.hist(angles, bins=20)
    plt.rc('text', usetex=True)
    plt.title(f"Histogram of the angles between " + r"$V_{th}$ and $V$" + 
              f" (T = {Thresholding_string})", fontsize = 16)
    plt.xlabel("Degrees", fontsize = 16)
    plt.axis([0, 180, 0, 90000])
    plt.savefig(f'{folder}Angle_Plots/{timestamp}_NA-{NODES}_AngleHist_{Thresholding_string.replace(".", "_")}.png')
    plt.rc('text', usetex=False)


    return VthTPxV

def calc_identity_err(VthTPxV):

    print("Raw Vth' * V")
    print(VthTPxV)

    err_austin = calc_identity_err_austin_method(VthTPxV)
    err_cm = calc_identity_err_cuthill_mckee_method(VthTPxV)

    return err_austin, err_cm

def calc_identity_err_cuthill_mckee_method(VthTPxV):
    sym = False
    VthTPxV_csr = csr_matrix(VthTPxV)
    perm = reverse_cuthill_mckee(VthTPxV_csr, sym)[::-1]
    
    for i in range(len(perm)):
        VthTPxV_csr[:,i] = VthTPxV_csr[perm,i]
    # for i in range(len(perm)):
    #     VthTPxV_csr[i,:] = VthTPxV_csr[i,perm]

    VthTPxV_perm = VthTPxV_csr.toarray()
    

    I_error = find_identity_error(VthTPxV_perm)
    print("\nError off identity matrix (using Cuthill Mckee method)", I_error)
    print("\nReordered Vth' * V (using Cuthill Mckee method)")
    print(VthTPxV_perm)

    return I_error

def calc_identity_err_austin_method(VthTPxV):

    # Reorder rows to identity using max value in each row
    for i in range(len(VthTPxV)):

        # Skip if value on diagonal is already approx 1
        if abs(1 - VthTPxV[i][i]) < 0.0001:
            continue
        
        # find max value in the row
        max_val, max_index = VthTPxV[0][i], 0
        for j in range(i, len(VthTPxV)):
            if abs(VthTPxV[j][i]) > abs(max_val):
                max_val = VthTPxV[j][i]
                max_index = j
        
        # place row in correct place
        tmp = numpy.copy(VthTPxV[max_index])
        VthTPxV[max_index] = numpy.copy(VthTPxV[i])
        VthTPxV[i] = tmp

    I_error = find_identity_error(VthTPxV)
    print("\nError off identity matrix (using austin's method)", I_error)
    print("\nReordered Vth' * V (using austin's method)")
    print(VthTPxV)

    return I_error

def find_identity_error(VthTPxV):
    # Find Error 
    I_error_sum1 = numpy.linalg.norm(numpy.abs(VthTPxV) - numpy.identity(NODES))
    I_error = I_error_sum1 / (NODES * NODES)
    return I_error



def compare_treshold_to_orig(val_orig, vec_orig, M_tresh, name):

    val_tresh, vec_tresh = get_evd_of_graph(M_tresh)
    get_tresholding_error(vec_orig, vec_tresh, name)
    return

    # depreciated
    x = numpy.linspace(0, NODES-1, NODES)


    fig, ((ax1, ax2)) = plt.subplots(2,1)
    ax1.set_title("Compare Eigenvalues and Eigenvectors" +
                 f" after Treshold = {name}")
    ax1.plot(x, val_orig, label = 'Original Eigenvalues')
    ax1.plot(x, val_tresh, label = 'Tresholded Eigenvalues')
    ax1.legend()
    ax1.set_ylabel("Eigenvalue")

    
    cos_sim = get_cos_sim(vec_orig, vec_tresh, name)
    ax2.plot(x, cos_sim)
    ax2.set_ylabel("Cosine Similarity")
    ax2.set_xlabel("Eigenvalue sorted smallest to biggest")
    ax2.set_ylim(-1.1, 1.1)
    fig.savefig(f'compare_tresh_{name}.png')


def print_graphs():

    global A, M, gft_graph_0_00005, gft_graph_0_005, gft_graph_0_05, gft_graph_0_1, gft_graph_0_2, gft_graph_0_3, gft_graph_1
    
    try:
        print("A")
        print(nx.to_numpy_array(A))
        print("\nM")
        print(nx.to_numpy_array(M))
        print("\ngft_graph_0_00005")                
        print(nx.to_numpy_array(gft_graph_0_00005))
        print("\ngft_graph_0_005")
        print(nx.to_numpy_array(gft_graph_0_005))
        print("\ngft_graph_0_05")               
        print(nx.to_numpy_array(gft_graph_0_05)) 
        print("\ngft_graph_0_1")                    
        print(nx.to_numpy_array(gft_graph_0_1))  
        print("\ngft_graph_0_2")                
        print(nx.to_numpy_array(gft_graph_0_2)) 
        print("\ngft_graph_0_3")                
        print(nx.to_numpy_array(gft_graph_0_3))    
        print("\ngft_graph_1")              
        print(nx.to_numpy_array(gft_graph_1))
    except:
        pass   
    
    quit()







###############################################
##                 Start Code                ##
###############################################


(make_dir("GSP", True, True))
(make_dir(folder + "Angle_Plots"))
(make_dir(folder + "Max_Val_Plots"))
redirect_to_file("log") 


# generate connected graph
A = gen_conn_graph()
adj_mtx = nx.to_numpy_array(A)
# print(adj_mtx)
# graph  = nx.from_numpy_array(adj_mtx + 1e-16)

# draw graph
plt.figure()
# nx.draw(A, pos=nx.circular_layout(A))
plt.title(f"Original A Graph (N(A) = {NODES})", fontsize = 16)
draw_graph(A, "A", SHOW_GRAPHS)
plt.savefig(f'{folder}{timestamp}_NA-{NODES}_A.png')
# plt.show()

#### get M ####

M, GFT, h, GFT_INV = get_gft_graph(A)
M = threshold(M, 1e-11) # remove false edges
print("\nA")
print(M)

A_reconstruct =nx.from_numpy_matrix(numpy.matmul(numpy.matmul(GFT_INV, h), GFT))
A_reconstruct = threshold(A_reconstruct, 1e-11) # remove false edges
print("\n\n")

plt.figure()
draw_graph(A_reconstruct, "A recon", SHOW_GRAPHS)
plt.title("Reconstructed A Graph", fontsize = 16)
plt.savefig(f'{folder}{timestamp}_A_recon.png')
# plt.show()

plt.figure()
plt.title(f"Original M Graph (N(A) = {NODES})", fontsize = 16)
draw_graph(M, "M", SHOW_GRAPHS)
plt.savefig(f'{folder}{timestamp}_NA-{NODES}_M.png')
# plt.show()

#### Thresholding ####

# get M evd
eval_orig, evec_orig = get_evd_of_graph(M)
compare_treshold_to_orig(eval_orig, evec_orig, M, "0")

threshold_factor = 0.00005
print(f"\nTresholding: {threshold_factor}")
gft_graph_0_00005 = threshold(M, threshold_factor)
plt.figure()
plt.title(f"M After Thresholding @ T = {threshold_factor}", fontsize = 16)
draw_graph(gft_graph_0_00005, f"Thresholded @ {threshold_factor}", SHOW_GRAPHS)
plt.savefig(f'{folder}{timestamp}_NA-{NODES}_M_threshold_0_00005.png')
# plt.show()
compare_treshold_to_orig(eval_orig, evec_orig, gft_graph_0_00005, f"{threshold_factor}")

threshold_factor = 0.005
print(f"\nTresholding: {threshold_factor}")
gft_graph_0_005 = threshold(M, threshold_factor)
plt.figure()
plt.title(f"M After Thresholding @ T = {threshold_factor}", fontsize = 16)
draw_graph(gft_graph_0_005, f"Thresholded @ {threshold_factor}", True)
plt.savefig(f'{folder}{timestamp}_NA-{NODES}_M_threshold_0_005.png')
# plt.show()
compare_treshold_to_orig(eval_orig, evec_orig, gft_graph_0_005, f"{threshold_factor}")


threshold_factor = 0.05
print(f"\nTresholding: {threshold_factor}")
gft_graph_0_05 = threshold(M, threshold_factor)
plt.figure()
plt.title(f"M After Thresholding @ T = {threshold_factor}", fontsize = 16)
draw_graph(gft_graph_0_05, f"Thresholded @ {threshold_factor}", True)
plt.savefig(f'{folder}{timestamp}_NA-{NODES}_M_threshold_0_05.png')
# plt.show()
compare_treshold_to_orig(eval_orig, evec_orig, gft_graph_0_05, f"{threshold_factor}")

threshold_factor = 0.1
print(f"\nTresholding: {threshold_factor}")
gft_graph_0_1 = threshold(M, threshold_factor)
plt.figure()
plt.title(f"M After Thresholding @ T = {threshold_factor}", fontsize = 16)
draw_graph(gft_graph_0_1, f"Thresholded @ {threshold_factor}", True)
plt.savefig(f'{folder}{timestamp}_NA-{NODES}_M_threshold_0_1.png')
# plt.show()
compare_treshold_to_orig(eval_orig, evec_orig, gft_graph_0_1, f"{threshold_factor}")

threshold_factor = 0.2
print(f"\nTresholding: {threshold_factor}")
gft_graph_0_2 = threshold(M, threshold_factor)
plt.figure()
plt.title(f"M After Thresholding @ T = {threshold_factor}", fontsize = 16)
draw_graph(gft_graph_0_2, f"Thresholded @ {threshold_factor}", True)
plt.savefig(f'{folder}{timestamp}_NA-{NODES}_M_threshold_0_2.png')
# plt.show()
compare_treshold_to_orig(eval_orig, evec_orig, gft_graph_0_2, f"{threshold_factor}")

threshold_factor = 0.3
print(f"\nTresholding: {threshold_factor}")
gft_graph_0_3 = threshold(M, threshold_factor)
plt.figure()
plt.title(f"M After Thresholding @ T = {threshold_factor}", fontsize = 16)
draw_graph(gft_graph_0_3, f"Thresholded @ {threshold_factor}", True)
plt.savefig(f'{folder}{timestamp}_NA-{NODES}_M_threshold_0_3.png')
# plt.show()
compare_treshold_to_orig(eval_orig, evec_orig, gft_graph_0_3, f"{threshold_factor}")


threshold_factor = 1
print(f"\nTresholding: {threshold_factor}")
gft_graph_1 = threshold(M, threshold_factor)
plt.figure()
plt.title(f"M After Thresholding @ T = {threshold_factor:1.5f}")
draw_graph(gft_graph_1, f"Thresholded @ {threshold_factor}", True)
plt.savefig(f'{folder}{timestamp}_NA-{NODES}_M_threshold_1.png')
# plt.show()
compare_treshold_to_orig(eval_orig, evec_orig, gft_graph_1, f"{threshold_factor}")



print_graphs()



