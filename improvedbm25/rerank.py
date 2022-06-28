import pymongo


myclient =pymongo.MongoClient("mongodb://localhost:27017/")

mytop=myclient['cs2019_top_1000']
myrerank=myclient['cs2019_rerank']
myrerank1000=myclient['cs2019_rerank_1000']
mydb = myclient["pubmed"]
mypapers=mydb["freqwords3"]#pubmed中文献信息表


def rerank(k1,k2,b1,b2,gg):
    for j in range(40):
        PMID_list=[]
        temp1 = "cs2019_top_1000_" + str(j + 1)
        mytop_100 = mytop[temp1]
        temp="cs2019_rerank_" + str(j + 1)
        rerank=myrerank[temp]
        i=0
        # for x in mytop_100.find({}, {'PMID', 'related', 'idf_para', 'cmk_len', 'cmk_freq', 'gx', 'bm25_score'}):
        #         PMID_list.append(x['PMID'])
        # print(PMID_list)
        # PMID_dict = {}
        # for y in mypapers.find({'PMID': {'$in': PMID_list}}, {'PMID', 'abstract'}):
        #     PMID_dict[y['PMID']] = y['abstract']
        for x in mytop_100.find({}, {'PMID', 'related', 'idf_para', 'cmk_len', 'cmk_freq', 'gx', 'bm25_score','ab_len','ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'}):
            idf_para_sum = 0
            idf_para = x['idf_para']
            cmk_len = x['cmk_len']
            cmk_freq = x['cmk_freq']
            gx = gg*float(x['gx'])
            ab_len = x['ab_len']
            for dic in idf_para:
                for key in dic:
                    idf_para_sum = idf_para_sum + (((k1 + 1) * float(key)) / (
                                (k1 * (b1 + (1 - b1) * (ab_len / 83))) + float(key))) * float(dic[key])
            cmk_socre = ((k2 + 1) * cmk_freq) / ((k2 * (b2 + (1 - b2) * (cmk_len / 13))) + cmk_freq)
            bm25_score = idf_para_sum + gx+cmk_socre

            mydict = {"PMID": x['PMID'], "related": x['related'],'ab_score':idf_para_sum, "idf_para": x['idf_para'],"ab_len":ab_len,
                      "cmk_len": x['cmk_len'], "cmk_freq": x['cmk_freq'],
                      "gx": x['gx'], "bm25_score": bm25_score,'ChemicalNameList':x['ChemicalNameList'], 'MeshHeadingNameList':x['MeshHeadingNameList'], 'KeywordsList':x['KeywordsList']}
            ss = rerank.insert_one(mydict)
            i=i+1
            print(str(ss)+'---------'+str(i))



if __name__ == '__main__':
        # rerank(3.21912, 98.66573,0.89633, 0.98359, 4.31791)
        # cal_p10()
        rerank(3.51495, 91.27053,0.84409, 1.00000, 4.03992)