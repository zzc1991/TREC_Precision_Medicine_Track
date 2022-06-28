import pymongo
import re
from math import log


myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["bert"]
mywords = mydb["freqwords1"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2019"]#pubmed中的主题词相关文献列表


mytopicdb=myclient["cs2019_bert"]
mydata=mytopicdb["cs2019_score_38"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2019_score_38_related"]#聚类后对应与主题相关联的文献


def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1=1.2
    k2=1.2
    b1=0.75
    b2=0.75
    idf_peutzjeghers = log((29138919 - 1868 + 0.5) / (1868 + 0.5), 10)
    idf_syndrome = log((29138919 - 891390 + 0.5) / (891390 + 0.5), 10)
    idf_stk11 = log((29138919 - 310 + 0.5) / (310 + 0.5), 10)


    idf_ele_1 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_2 = log((13670358 - 858+ 0.5) / (858 + 0.5), 10)
    idf_ele_3 = log((13670358 - 56231 + 0.5) / (56231 + 0.5), 10)
    idf_ele_4 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_5 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_6 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_7 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)


    idf_eleM_1 = log((25389659 - 1848 + 0.5) / (1848 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 0 + 0.5) / (0+ 0.5), 10)
    idf_eleM_3 = log((25389659 - 56231 + 0.5) / (56231 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 8002162 + 0.5) / (8002162 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 1862883 + 0.5) / (1862883 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 248 + 0.5) / (248 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_11 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 112 + 0.5) / (112 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 63 + 0.5) / (63 + 0.5), 10)
    idf_eleK_3 = log((5435471 - 0 + 0.5) / (0 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        gx = 0
        gx1 = 0
        gx2 = 0
        gx3 = 0
        gx4=0
        len_freq=0
        peutzjeghers_score=0
        syndrome_score = 0
        stk11_score = 0
        cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
        ChemicalNameList = x['ChemicalNameList']
        MeshHeadingNameList = x['MeshHeadingNameList']
        KeywordsList = x['KeywordsList']
        wordfreq = x['wordfreq']
        peutzjeghers = [True for x in wordfreq.items() if 'peutz-jeghers' in x]
        syndrome = [True for x in wordfreq.items() if 'syndrome' in x]

        # ---------------摘要统计-------------------#


        for key in wordfreq:
            len_freq = len_freq + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'peutz-jeghers' in key1:
                peutzjeghers_score = peutzjeghers_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'syndrome' in key1:
                syndrome_score = syndrome_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'stk11' in key1:
                stk11_score = stk11_score + wordfreq[key]



        bm25_peutzjeghers_score = (((k1+1)*peutzjeghers_score)/((k1*(b1+(1-b1)*(len_freq/85)))+peutzjeghers_score))
        bm25_syndrome_score = (((k1 + 1) * syndrome_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + syndrome_score))
        bm25_stk11_score = (((k1 + 1) * stk11_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + stk11_score))

        bm25_ab_score =idf_peutzjeghers*bm25_peutzjeghers_score+idf_syndrome*bm25_syndrome_score+idf_stk11*bm25_stk11_score

        idf_para=[{str(peutzjeghers_score):idf_peutzjeghers},{str(syndrome_score):idf_syndrome},{str(stk11_score):idf_stk11}]

        # ---------------共现分析摘要-------------------#
        if len(peutzjeghers) != 0 and peutzjeghers[0] and len(syndrome) != 0 and syndrome[0]:
            for key in wordfreq:
                key = cop.sub('', key)
                if 'stk11' in key:
                    gx = idf_stk11

    # ---------------共现分析化学-------------------#
        if len(peutzjeghers) != 0 and peutzjeghers[0] and len(syndrome) != 0 and syndrome[0]:
            for ele in ChemicalNameList:
                if 'STK11' in ele['NameOfSubstance']:
                    gx = idf_stk11
                    break

    # ---------------共现分析关键字-------------------#
        if len(peutzjeghers) != 0 and peutzjeghers[0] and len(syndrome) != 0 and syndrome[0]:
            for eleK in KeywordsList:
                if 'stk11' in str(eleK).lower():
                    gx = idf_stk11
                    break

     # ---------------共现分析医学主题词-------------------#
        if len(peutzjeghers) != 0 and peutzjeghers[0] and len(syndrome) != 0 and syndrome[0]:
            for eleM in MeshHeadingNameList:
                if 'STK11' in eleM['MeshHeadingName']:
                    gx = idf_stk11
                    break


        for ele in ChemicalNameList:
            if 'Peutz-Jeghers Syndrome' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_1
                break
        for ele in ChemicalNameList:
            if 'STK11 protein, human' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_2
                break
        for ele in ChemicalNameList:
            if 'Protein-Serine-Threonine Kinases' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_3
                break





        for eleM in MeshHeadingNameList:
            if 'Peutz-Jeghers Syndrome' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_1
                break
        for eleM in MeshHeadingNameList:
            if 'STK11 protein, human' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_2
                break
        for eleM in MeshHeadingNameList:
            if 'Protein-Serine-Threonine Kinases' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_3
                break

        for eleM in MeshHeadingNameList:
            if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_5
                break
        for eleM in MeshHeadingNameList:
            if 'Male' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_6
                break
        for eleM in MeshHeadingNameList:
            if 'Child' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_7
                break



        for eleK in KeywordsList:
            if 'peutz-jeghers syndrome' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_1
                break
        for eleK in KeywordsList:
            if 'stk11' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_2
                break


        total_gx=gx1+gx2+gx3+gx+gx4
        cmk_len=len(ChemicalNameList) + len(MeshHeadingNameList) + len(KeywordsList)
        bm25_cmk_len=ss1 + ss2 + ss4
        bm25_cmk_score = (((k2 + 1) * bm25_cmk_len) / ((k2 * (b2 + (1 - b2) * (cmk_len / 13))) + bm25_cmk_len))
        bm25_score=bm25_ab_score+bm25_cmk_score+total_gx
        if(bm25_score>yuzhi):
            mydict = {"PMID": x['PMID'],"ab_score":bm25_ab_score,"idf_para":idf_para,
                      "cmk_len":cmk_len,"cmk_freq":bm25_cmk_len,"bm25_cmk_score":bm25_cmk_score,"gx":total_gx,"bm25_score":bm25_score,
                       "ChemicalNameList":x['ChemicalNameList'],"MeshHeadingNameList":x['MeshHeadingNameList'],"KeywordsList":x['KeywordsList']}
            y = mydata.insert_one(mydict)
            k=k+1
            print(str(y) + '---------' + str(k))

def count(mysort,mycount,topic):
    for x in mysort.find({}, {'PMID', 'ab_score','idf_para', 'cmk_len', 'cmk_freq', 'bm25_cmk_score','gx','bm25_score',
                              'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'}):
        kk = 0
        for y in mytopic.find({"topic": topic}, {'PMID', 'relate'}):
            if x['PMID'] == y['PMID']:
                mydict = {"PMID": x['PMID'], "related": y['relate'], "ab_score":x["ab_score"],"idf_para":x['idf_para'],
                          "cmk_len": x['cmk_len'], "cmk_freq": x['cmk_freq'],'bm25_cmk_score':x['bm25_cmk_score'],'gx':x['gx'],
                          "bm25_score": x['bm25_score'],
                          "ChemicalNameList": x['ChemicalNameList'], "MeshHeadingNameList": x['MeshHeadingNameList'],
                          "KeywordsList": x['KeywordsList']}
                ss = mycount.insert_one(mydict)
                print(ss)
                kk = kk + 1
        if (kk == 0):
            mydict = {"PMID": x['PMID'], "related": -1, "ab_score": x["ab_score"], "idf_para": x['idf_para'],
                      "cmk_len": x['cmk_len'], "cmk_freq": x['cmk_freq'], 'bm25_cmk_score': x['bm25_cmk_score'],
                      'gx': x['gx'],
                      "bm25_score": x['bm25_score'],
                      "ChemicalNameList": x['ChemicalNameList'], "MeshHeadingNameList": x['MeshHeadingNameList'],
                      "KeywordsList": x['KeywordsList']}
            ss = mycount.insert_one(mydict)
            print(ss)

if __name__ == '__main__':
    sortsecond(mywords,mydata,4)
    count(mydata,mycount,"38")



