import pymongo
import pyterrier.utils as ptu
import pyterrier.io as pto
import pandas as pd
import numpy as np

myclient = pymongo.MongoClient("mongodb://localhost:27017/")
qrels = pto.read_qrels('PM_result/qrels-treceval-abstracts.2019.txt')
mytopicdb = myclient["cs2019_nomalize"]
res = []

number=40
for i in range(number):
    conut=1
    ss='result_add_'+str(i+1)
    top_1000=mytopicdb[ss]
    for doc in top_1000.find().sort('score', -1).limit(500):
        res.append([i+1, doc['PMID'], conut, doc['com_score']])
        conut += + 1

res = pd.DataFrame(res, columns=['qid', 'docno', 'rank', 'score'])
pto.write_results(res, 'result')
res=pto.read_results('result')

result=ptu.Utils.evaluate(res, qrels, metrics = ['map','ndcg','P','Rprec','recall',])
print(result)


# score_dict={}
# for i in range(number):
#     score=[]
#     ss='cs2019_'+str(i+1)+'_top1000_s'
#     top_1000=mytopicdb[ss]
#     for doc in top_1000.find().sort('score', -1):
#         score.append(float(doc['score']))
#     score_dict[i+1]=score
#
# for item in score_dict:
#     print(item,np.mean(score_dict[item]),np.std(score_dict[item]),np.median(score_dict[item]))
