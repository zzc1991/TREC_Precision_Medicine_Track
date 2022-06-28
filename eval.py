import csv
from collections import Counter

import pymongo
import pyterrier.utils as ptu
import pyterrier.io as pto
import pandas as pd
import numpy as np

myclient = pymongo.MongoClient("mongodb://localhost:27017/")
qrels = pto.read_qrels('PM_result/qrels-treceval-abstracts.2019.txt')
mytopicdb = myclient["cs2019_nomalize"]
pubmed = myclient["pubmed"]
topic=pubmed['topics2019']
res = []
all_static = []
not_ex = 0
number = 40

f1=open('huatu.csv', 'w', encoding='utf-8', newline='')
# 2. 基于文件对象构建 csv写入对象
csv_writer = csv.writer(f1)
r_num=[]
for i in range(40):
    num=0
    for x in topic.find({'topic':str(i+1),'relate':{'$in':['1','2']}}):
        num+=1
    r_num.append(num)
print(r_num)
for i in range(number):
    static = []
    conut = 1
    # ss='cs2019_'+str(i+1)+'_top1000'
    ss = 'result' + str(i + 1) + '_all'
    top_1000 = mytopicdb[ss]
    for doc in top_1000.find().sort('score',-1).limit(1000):
        if conut <= r_num[i]:
            if int(doc['related']) == 1 or int(doc['related']) == 2:
                static.append(1)
            elif int(doc['related']) == 0:
                static.append(0)
            elif int(doc['related']) == -1:
                static.append(-1)
        res.append([i + 1, doc['PMID'], conut, doc['score']])
        conut += + 1
    all_static.append(Counter(static))


p10_list=[]
for i,item in enumerate(all_static):
    if 1 in item.keys():
        p10_list.append(item[1]/r_num[i])
    else:
        p10_list.append(0)

print(p10_list)

print(all_static)
print(not_ex)
res = pd.DataFrame(res, columns=['qid', 'docno', 'rank', 'score'])
pto.write_results(res, 'result')
res = pto.read_results('result')

result = ptu.Utils.evaluate(res, qrels, metrics=['map', 'P', 'Rprec', 'recall', 'ndcg'])
print(result)

