import pymongo
from math import log
myclient =pymongo.MongoClient("mongodb://localhost:27017/")


mytop=myclient['cs2019_rerank']
myres=myclient['cs2019_result']
def top():
    k=0
    for j in range(40):
        temp0="cs2019_rerank_"+str(j+1)
        mytopicdb=mytop[temp0]
        temp="cs2019_result_"+str(j+1)
        mydata = myres[temp]
        for x in mytopicdb.find({},
                             {'PMID', 'ab_score', 'bm25_cmk_score', 'related', 'idf_para', 'cmk_len', 'cmk_freq', 'gx',
                              'bm25_score', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'}).sort( [("bm25_score", -1)]):
                mydict = {"PMID": x['PMID'], "related": x['related'], "idf_para": x['idf_para'],
                          "bm25_score": x['bm25_score'],
                          "gx": x['gx'],
                          'ChemicalNameList': x['ChemicalNameList'], "cmk_len": x['cmk_len'], "cmk_freq": x['cmk_freq'],
                          'MeshHeadingNameList': x['MeshHeadingNameList'], 'KeywordsList': x['KeywordsList']}
                ss = mydata.insert_one(mydict)
                k = k + 1
                print(str(ss) + '-------------' + str(k))


if __name__ == '__main__':
     #statis(mywords,15,myfreq)
    top()
    pass
