import pymongo
import  pandas as pd
import math
import numpy as np
import scipy.special as sc_special
from math import log

myclient =pymongo.MongoClient("mongodb://localhost:27017/")

mytop=myclient['cs2019_top_1000']

def cuckoo_search(n, m, fit_func, lower_boundary, upper_boundary, iter_num=100, pa=0.25, beta=1.5, step_size=0.1,ite_numer=0):
    """
    Cuckoo search function
    ---------------------------------------------------
    Input parameters:
        n: Number of nests
        m: Number of dimensions
        fit_func: User defined fitness evaluative function
        lower_boundary: Lower bounary (example: lower_boundary = (-2, -2, -2))
        upper_boundary: Upper boundary (example: upper_boundary = (2, 2, 2))
        iter_num: Number of iterations (default: 100)
        pa: Possibility that hosts find cuckoos' eggs (default: 0.25)
        beta: Power law index (note: 1 < beta < 2) (default: 1.5)
        step_size:  Step size scaling factor related to the problem's scale (default: 0.1)
    Output:
        The best solution and its value
    """
    # get initial nests' locations

    nests = generate_nests(n, m, lower_boundary, upper_boundary)
    fitness = calc_fitness(fit_func, nests)

    # get the best nest and record it
    best_nest_index = np.argmax(fitness)
    best_fitness = fitness[best_nest_index]
    best_nest = nests[best_nest_index].copy()

    for _ in range(iter_num):
        ite_numer = ite_numer + 1
        print('iter==========' + str(ite_numer))
        nests = update_nests(fit_func, lower_boundary, upper_boundary, nests, best_nest, fitness, step_size)
        nests = abandon_nests(nests, lower_boundary, upper_boundary, pa)
        fitness = calc_fitness(fit_func, nests)

        max_nest_index = np.argmax(fitness)
        max_fitness = fitness[max_nest_index]
        max_nest = nests[max_nest_index]

        if (max_fitness > best_fitness):
            best_nest = max_nest.copy()
            best_fitness = max_fitness
        iter_list.append(best_fitness)
        print(best_fitness)
    return (best_nest, best_fitness)


def generate_nests(n, m, lower_boundary, upper_boundary):
    """
    Generate the nests' locations
    ---------------------------------------------------
    Input parameters:
        n: Number of nests
        m: Number of dimensions
        lower_boundary: Lower bounary (example: lower_boundary = (-2, -2, -2))
        upper_boundary: Upper boundary (example: upper_boundary = (2, 2, 2))
    Output:
        generated nests' locations
    """
    lower_boundary = np.array(lower_boundary)
    upper_boundary = np.array(upper_boundary)
    nests = np.empty((n, m))

    for each_nest in range(n):
        nests[each_nest] = lower_boundary + np.array([np.random.rand() for _ in range(m)]) * (
                    upper_boundary - lower_boundary)

    return nests


def update_nests(fit_func, lower_boundary, upper_boundary, nests, best_nest, fitness, step_coefficient):
    """
    This function is to get new nests' locations and use new better one to replace the old nest
    ---------------------------------------------------
    Input parameters:
        fit_func: User defined fitness evaluative function
        lower_boundary: Lower bounary (example: lower_boundary = (-2, -2, -2))
        upper_boundary: Upper boundary (example: upper_boundary = (2, 2, 2))
        nests: Old nests' locations
        best_nest: Nest with best fitness
        fitness: Every nest's fitness
        step_coefficient:  Step size scaling factor related to the problem's scale (default: 0.1)
    Output:
        Updated nests' locations
    """
    lower_boundary = np.array(lower_boundary)
    upper_boundary = np.array(upper_boundary)
    n, m = nests.shape
    # generate steps using levy flight
    steps = levy_flight(n, m, 1.5)
    new_nests = nests.copy()

    for each_nest in range(n):
        # coefficient 0.01 is to avoid levy flights becoming too aggresive
        # and (nest[each_nest] - best_nest) could let the best nest be remained
        step_size = step_coefficient * steps[each_nest] * (nests[each_nest] - best_nest)
        step_direction = np.random.rand(m)
        new_nests[each_nest] += step_size * step_direction
        # apply boundary condtions
        new_nests[each_nest][new_nests[each_nest] < lower_boundary] = lower_boundary[
            new_nests[each_nest] < lower_boundary]
        new_nests[each_nest][new_nests[each_nest] > upper_boundary] = upper_boundary[
            new_nests[each_nest] > upper_boundary]

    new_fitness = calc_fitness(fit_func, new_nests)
    nests[new_fitness > fitness] = new_nests[new_fitness > fitness]

    return nests


def abandon_nests(nests, lower_boundary, upper_boundary, pa):
    """
    Some cuckoos' eggs are found by hosts, and are abandoned.So cuckoos need to find new nests.
    ---------------------------------------------------
    Input parameters:
        nests: Current nests' locations
        lower_boundary: Lower bounary (example: lower_boundary = (-2, -2, -2))
        upper_boundary: Upper boundary (example: upper_boundary = (2, 2, 2))
        pa: Possibility that hosts find cuckoos' eggs
    Output:
        Updated nests' locations
    """
    lower_boundary = np.array(lower_boundary)
    upper_boundary = np.array(upper_boundary)
    n, m = nests.shape
    for each_nest in range(n):
        if (np.random.rand() < pa):
            step_size = np.random.rand() * (nests[np.random.randint(0, n)] - nests[np.random.randint(0, n)])
            nests[each_nest] += step_size
            # apply boundary condtions
            nests[each_nest][nests[each_nest] < lower_boundary] = lower_boundary[nests[each_nest] < lower_boundary]
            nests[each_nest][nests[each_nest] > upper_boundary] = upper_boundary[nests[each_nest] > upper_boundary]

    return nests


def levy_flight(n, m, beta):
    """
    This function implements Levy's flight.
    ---------------------------------------------------
    Input parameters:
        n: Number of steps
        m: Number of dimensions
        beta: Power law index (note: 1 < beta < 2)
    Output:
        'n' levy steps in 'm' dimension
    """
    sigma_u = (sc_special.gamma(1 + beta) * np.sin(np.pi * beta / 2) / (
                sc_special.gamma((1 + beta) / 2) * beta * (2 ** ((beta - 1) / 2)))) ** (1 / beta)
    sigma_v = 1

    u = np.random.normal(0, sigma_u, (n, m))
    v = np.random.normal(0, sigma_v, (n, m))

    steps = u / ((np.abs(v)) ** (1 / beta))

    return steps


def calc_fitness(fit_func, nests):
    """
    calculate each nest's fitness
    ---------------------------------------------------
    Input parameters:
        fit_func: User defined fitness evaluative function
        nests:  Nests' locations
    Output:
        Every nest's fitness
    """
    n, m = nests.shape
    fitness = np.empty(n)

    for each_nest in range(n):
        fitness[each_nest] = fit_func(nests[each_nest])

    return fitness
def IDCG(n):
    idcg=0
    for i in range(n):
        idcg+=1/math.log(i+2,2)
    return idcg

def nDCG(ranked_list,len):
    dcg=0
    idcg=IDCG(len)
    for i in range(len):
        if ranked_list[i] not in [1,2]:
            continue
        rank=i+1
        dcg+=1/math.log(rank+1,2)
    return  dcg/idcg if idcg!=0 else 0
def caculate(ranked_list):
    sum=0
    total=0
    for j in ranked_list[0:1000, 1]:
        if j in [1, 2]:
            sum = sum + 1
        if j in [1, 2,0,-1]:
            total = total + 1
        if total==10:
            break
    return sum/total
if __name__ == '__main__':
    iter_list=[]
    topic_mat_list=[]
    total_list=[]
    for j in range(40):
        score_mat = np.zeros(shape=(1, 2))
        temp1 = "cs2019_top_1000_" + str(j + 1)
        mytop_100 = mytop[temp1]
        list=[]
        df = pd.DataFrame(columns=['PMID', 'related', 'idf_para', 'cmk_len', 'cmk_freq', 'gx','bm25_score','ab_len'])
        for x in mytop_100.find({}, {'PMID', 'related', 'idf_para', 'cmk_len', 'cmk_freq', 'gx', 'bm25_score','ab_len'}).sort([("bm25_score", -1)]).limit(1000):
            list.append({'PMID':x['PMID'], 'related':x['related'], 'idf_para':x['idf_para'], 'cmk_len':x['cmk_len'], 'cmk_freq':x['cmk_freq'], 'gx':x['gx'],'bm25_score':x['bm25_score'],'ab_len':x['ab_len']})
            score_mat=np.row_stack((score_mat, [float(x['bm25_score']),float(x['related'])]))
        total_list.append(list)
        score_mat = np.delete(score_mat, 0, axis=0)
        topic_mat_list.append(score_mat)
        # print(score_mat)
    ss=total_list[0][0]
    print(ss)
    print(ss['idf_para'])
        #     score_mat=np.row_stack((score_mat, [float(x['ab_score']),float(x['Chemical_score']),float(x['Mesh_score']),float(x['Key_score']),float(x['gx']),0,x['related'],float(x['length'])]))
        # score_mat=np.delete(score_mat,0,axis=0)
        # topic_mat_list.append(score_mat)

    # score_mat = np.zeros(shape=(1, 7))
    # for x in mycount.find({}, {'PMID','ab_score','Chemical_score', 'Mesh_score', 'Key_score','gx','related'}):
    #     score_mat=np.row_stack((score_mat, [float(x['ab_score']),float(x['Chemical_score']),float(x['Mesh_score']),float(x['Key_score']),float(x['gx']),0,x['related']]))
    # score_mat=np.delete(score_mat,0,axis=0)
    #
    #
    # print( score_mat.shape[0])
    # for i in  range(score_mat.shape[0]):
    #       score_mat[i][5]=0.00778*float(score_mat[i][0])+0.01959*float(score_mat[i][1])+0.00131*float(score_mat[i][2])+0.06169*float(score_mat[i][3])+2.68604*float(score_mat[i][4])
    #        # score_mat[i][5] = -0.81476 * float(score_mat[i][0]) + -0.20987 * float(score_mat[i][1]) + -0.08953 * float(score_mat[i][2]) + -1.59468 * float(score_mat[i][3]) + 1.97928 * float(score_mat[i][4])
    # score_order = np.lexsort(score_mat.T[:6, :])
    # score_rev = score_order[::-1]
    # score_result = score_mat[score_rev, :]
    # ng=nDCG(score_result[:,6],len(score_result))
    # print(score_result[:,6])
    # for i in score_result[:,6]:
    #     print(i)
    #
    #
    def fit_func(nest):
        sum_ng=0
        k1, k2,b1,b2,ss= nest
        for k in range(40):
            score_mat=topic_mat_list[k]
            data_list=total_list[k]
            for i in range(score_mat.shape[0]):
                idf_para_sum=0
                idf_para=data_list[i]['idf_para']
                cmk_len = data_list[i]['cmk_len']
                cmk_freq = data_list[i]['cmk_freq']
                gx =data_list[i]['gx']
                ab_len = data_list[i]['ab_len']
                for dic in idf_para:
                    for key in dic:
                         idf_para_sum=idf_para_sum+(((k1+1)*float(key))/((k1*(b1+(1-b1)*(ab_len/85)))+float(key)))*float(dic[key])
                cmk_socre=((k2+1)*cmk_freq)/((k2*(b2+(1-b2)*(cmk_len/13)))+cmk_freq)
                bm25_score=idf_para_sum+cmk_socre+ss*gx
                score_mat[i][0] =bm25_score
            score_order = np.lexsort(score_mat.T[:1, :])
            score_rev = score_order[::-1]
            score_result = score_mat[score_rev, :]
            sum=caculate(score_result)
            # sum_ng =sum_ng+nDCG(score_result[0:100, 1], 100)
            sum_ng = sum_ng+sum
        return sum_ng/40
    best_nest, best_fitness = cuckoo_search(40, 5, fit_func, [0,0,0,0,0], [100, 100,1,1,5], step_size = 0.4,iter_num=500)
    print('最大值为:%.5f, 在(%.5f, %.5f,%.5f, %.5f, %.5f)处取到!' % (best_fitness, best_nest[0],best_nest[1], best_nest[2],best_nest[3],best_nest[4]))
    print(iter_list)
