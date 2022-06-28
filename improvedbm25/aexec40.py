import pymongo
import re
from math import log


myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["bert"]
mywords = mydb["freqwords1"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2019"]#pubmed中的主题词相关文献列表


mytopicdb=myclient["cs2019_bert"]
mydata=mytopicdb["cs2019_score_40"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2019_score_40_related"]#聚类后对应与主题相关联的文献


def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1=1.2
    k2=1.2
    b1=0.75
    b2=0.75
    idf_malignant = log((29138919 - 358128 + 0.5) / (358128 + 0.5), 10)
    idf_hyperthermia = log((29138919 - 23780 + 0.5) / (23780 + 0.5), 10)
    idf_ryr1 = log((29138919 - 1391 + 0.5) / (1391 + 0.5), 10)


    idf_ele_1 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0+ 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_4 = log((13670358 - 5703 + 0.5) / (5703 + 0.5), 10)
    idf_ele_5 = log((13670358 - 259884 + 0.5) / (259884 + 0.5), 10)
    idf_ele_6 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_7 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)


    idf_eleM_1 = log((25389659 - 3559 + 0.5) / (3559 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 93956 + 0.5) / (93956+ 0.5), 10)
    idf_eleM_3 = log((25389659 - 131193 + 0.5) / (131193 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 5703+ 0.5) / (5703 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 8002162 + 0.5) / (8002162 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 4029038 + 0.5) / (4029038 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 259884 + 0.5) / (259884+ 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_11 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 157 + 0.5) / (157 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 74 + 0.5) / (74 + 0.5), 10)
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
        malignant_score=0
        hyperthermia_score = 0
        ryr1_score = 0
        cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
        ChemicalNameList = x['ChemicalNameList']
        MeshHeadingNameList = x['MeshHeadingNameList']
        KeywordsList = x['KeywordsList']
        wordfreq = x['wordfreq']
        malignant= [True for x in wordfreq.items() if 'malignant' in x]
        hyperthermia = [True for x in wordfreq.items() if 'hyperthermia' in x]

        # ---------------摘要统计-------------------#


        for key in wordfreq:
            len_freq = len_freq + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'malignant' in key1:
                malignant_score = malignant_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'hyperthermia' in key1:
                hyperthermia_score = hyperthermia_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'ryr1' in key1:
                ryr1_score = ryr1_score + wordfreq[key]



        bm25_malignant_score = (((k1+1)*malignant_score)/((k1*(b1+(1-b1)*(len_freq/85)))+malignant_score))
        bm25_hyperthermia_score = (((k1 + 1) * hyperthermia_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + hyperthermia_score))
        bm25_ryr1_score = (((k1 + 1) * ryr1_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + ryr1_score))

        bm25_ab_score =idf_malignant*bm25_malignant_score+idf_hyperthermia*bm25_hyperthermia_score+idf_ryr1*bm25_ryr1_score

        idf_para=[{str(malignant_score):idf_malignant},{str(hyperthermia_score):idf_hyperthermia},{str(ryr1_score):idf_ryr1}]

        # ---------------共现分析摘要-------------------#
        if len(malignant) != 0 and malignant[0] and len(hyperthermia) != 0 and hyperthermia[0]:
            for key in wordfreq:
                key = cop.sub('', key)
                if 'ryr1' in key:
                    gx = idf_ryr1

    # ---------------共现分析化学-------------------#
        if len(malignant) != 0 and malignant[0] and len(hyperthermia) != 0 and hyperthermia[0]:
            for ele in ChemicalNameList:
                if 'RYR1' in ele['NameOfSubstance']:
                    gx = idf_ryr1
                    break

    # ---------------共现分析关键字-------------------#
        if len(malignant) != 0 and malignant[0] and len(hyperthermia) != 0 and hyperthermia[0]:
            for eleK in KeywordsList:
                if 'ryr1' in str(eleK).lower():
                    gx = idf_ryr1
                    break

     # ---------------共现分析医学主题词-------------------#
        if len(malignant) != 0 and malignant[0] and len(hyperthermia) != 0 and hyperthermia[0]:
            for eleM in MeshHeadingNameList:
                if 'RYR1' in eleM['MeshHeadingName']:
                    gx = idf_ryr1
                    break


        for ele in ChemicalNameList:
            if 'Malignant Hyperthermia' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_1
                break
        for ele in ChemicalNameList:
            if 'Muscle Contraction' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_2
                break
        for ele in ChemicalNameList:
            if 'Muscle, Skeletal' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_3
                break
        for ele in ChemicalNameList:
            if 'Ryanodine Receptor Calcium Release Channel' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_4
                break
        for ele in ChemicalNameList:
            if 'Calcium' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_5
                break




        for eleM in MeshHeadingNameList:
            if 'Malignant Hyperthermia' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_1
                break
        for eleM in MeshHeadingNameList:
            if 'Muscle Contraction' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_2
                break
        for eleM in MeshHeadingNameList:
            if 'Muscle, Skeletal' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_3
                break
        for eleM in MeshHeadingNameList:
            if 'Ryanodine Receptor Calcium Release Channel' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_4
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
            if 'Middle Aged' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_7
                break
        for eleM in MeshHeadingNameList:
            if 'Calcium' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_8
                break
        for eleM in MeshHeadingNameList:
            if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_9
                break

        for eleK in KeywordsList:
            if 'malignant hyperthermia' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_1
                break
        for eleK in KeywordsList:
            if 'ryr1' in str(eleK).lower():
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
    sortsecond(mywords,mydata,9)
    count(mydata,mycount,"40")



