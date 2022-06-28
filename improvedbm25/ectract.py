import pymongo
from math import log
myclient =pymongo.MongoClient("mongodb://localhost:27017/")


mytop=myclient['cs2019_top_1000']
mydb = myclient["pubmed"]
mypapers=mydb["freqwords3"]#pubmed中文献信息表
def top(number):
    for j in range(40):
        temp0="cs2019_"+str(j+1)
        mytopicdb=myclient[temp0]
        temp="cs2019_score_"+str(j+1)+"_related"
        mydata = mytopicdb[temp]
        temp1="cs2019_top_1000_"+str(j+1)
        mytop_100=mytop[temp1]
        i=0
        k=0
        PMID_list=[]
        wordfreq_list=[]
        for x in mydata.find({}, {'PMID', 'related', 'idf_para', 'cmk_len', 'cmk_freq', 'gx', 'bm25_score'}).sort([("bm25_score",-1)]):
            if i<number:
                PMID_list.append(x['PMID'])
                # mydict = {"PMID": x['PMID'], "related": x['related'], "idf_para": x['idf_para'],
                #           "cmk_len": x['cmk_len'], "cmk_freq": x['cmk_freq'],
                #           "gx": x['gx'], "bm25_score": x['bm25_score'],}
                # ss = mytop_100.insert_one(mydict)
                # print(str(ss)+'-------------'+str(i))
                i=i+1
        print(PMID_list)
        dic = {}
        for y in mypapers.find({'PMID': {'$in': PMID_list}}, {'PMID', 'wordfreq'}):
            len_freq=0
            wordfreq = y['wordfreq']
            for key in wordfreq:
                len_freq = len_freq + wordfreq[key]
            dic[y['PMID']]=len_freq
        print(dic)
        for x in mydata.find({}, {'PMID','ab_score','bm25_cmk_score','related', 'idf_para', 'cmk_len', 'cmk_freq', 'gx', 'bm25_score','ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'}).sort([("bm25_score", -1)]):
            if x['PMID'] in dic:
                gx=log(float(x['gx'])+1,2)
                bm25_score=float(x['ab_score'])+float(x['bm25_cmk_score'])+gx
                mydict = {"PMID": x['PMID'], "related": x['related'], "idf_para": x['idf_para'],
                          "bm25_score": x['bm25_score'], "bm25log_score": bm25_score,
                          "gx": x['gx'],'log_gx':gx,"ab_len":dic[x['PMID']], 'ChemicalNameList':x['ChemicalNameList'],"cmk_len": x['cmk_len'], "cmk_freq": x['cmk_freq'], 'MeshHeadingNameList':x['MeshHeadingNameList'], 'KeywordsList':x['KeywordsList']}
                ss = mytop_100.insert_one(mydict)
                k = k + 1
                print(str(ss)+'-------------'+str(k))

if __name__ == '__main__':
     #statis(mywords,15,myfreq)
    top(1000)
    pass
