import pymongo
import re
from math import log


myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["pubmed"]
mywords = mydb["freqwords1"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2019"]#pubmed中的主题词相关文献列表


mytopicdb=myclient["cs2019_bert"]
mydata=mytopicdb["cs2019_score_2"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2019_score_2_related"]#聚类后对应与主题相关联的文献





def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1=1.2
    k2=1.2
    b1=0.75
    b2=0.75
    idf_cholangiocarcinoma = log((29138919 - 9173 + 0.5) / (9173 + 0.5), 10)
    idf_braf = log((29138919 - 8527 + 0.5) / (8527 + 0.5), 10)
    idf_v600e = log((29138919 - 1887 + 0.5) / (1887 + 0.5), 10)

    idf_ele_1 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_4 = log((13670358 - 7222 + 0.5) / (7222 + 0.5), 10)
    idf_ele_5 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_6 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)


    idf_eleM_1 = log((25389659 - 7716 + 0.5) / (7716 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 13396 + 0.5) / (13396 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 10185 + 0.5) / (10185 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 7221 + 0.5) / (7221 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 8002162 + 0.5) / (8002162 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 2842020 + 0.5) / (2842020 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 78561 + 0.5) / (78561 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 113103 + 0.5) / (113103 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 2052 + 0.5) / (2052 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 415 + 0.5) / (415 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        gx = 0
        gx1 = 0
        gx2 = 0
        gx3 = 0
        len_freq=0
        cholangiocarcinoma_score=0
        braf_score=0
        v600e_score=0

        cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
        ChemicalNameList = x['ChemicalNameList']
        MeshHeadingNameList = x['MeshHeadingNameList']
        KeywordsList = x['KeywordsList']
        wordfreq = x['wordfreq']
        cholangiocarcinoma = [True for x in wordfreq.items() if 'cholangiocarcinoma' in x]


        # ---------------摘要统计-------------------#


        for key in wordfreq:
            len_freq = len_freq + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'cholangiocarcinoma' in key1:
                cholangiocarcinoma_score = cholangiocarcinoma_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'braf' in key1:
                braf_score = braf_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'v600e' in key1:
                v600e_score = v600e_score + wordfreq[key]

        bm25_cholangiocarcinoma_score = (((k1+1)*cholangiocarcinoma_score)/((k1*(b1+(1-b1)*(len_freq/85)))+cholangiocarcinoma_score))
        bm25_braf_score =(((k1+1)*braf_score)/((k1*(b1+(1-b1)*(len_freq/85)))+braf_score))
        bm25_v600e_score = (((k1+1)*v600e_score)/((k1*(b1+(1-b1)*(len_freq/85)))+v600e_score))
        bm25_ab_score =idf_cholangiocarcinoma*bm25_cholangiocarcinoma_score+idf_braf*bm25_braf_score+idf_v600e*bm25_v600e_score

        idf_para=[{str(cholangiocarcinoma_score):idf_cholangiocarcinoma},{str(braf_score):idf_braf},{str(v600e_score):idf_v600e}]

        # ---------------共现分析摘要-------------------#
        if len(cholangiocarcinoma)!=0 and cholangiocarcinoma[0]:
            for key in wordfreq:
                key = cop.sub('', key)
                if 'braf' in key:
                    gx = idf_braf
                if 'v600e' in key:
                    gx1 = idf_v600e

    # ---------------共现分析化学-------------------#
        if len(cholangiocarcinoma) != 0 and cholangiocarcinoma[0]:
            for ele in ChemicalNameList:
                if 'V600E' in ele['NameOfSubstance']:
                    gx1 = idf_v600e
                    break
        if len(cholangiocarcinoma) != 0 and cholangiocarcinoma[0]:
            for ele in ChemicalNameList:
                if 'B-raf' in ele['NameOfSubstance']:
                    gx = idf_braf
                    break
    # ---------------共现分析关键字-------------------#
        if len(cholangiocarcinoma) != 0 and cholangiocarcinoma[0]:
            for eleK in KeywordsList:
                if 'v600e' in str(eleK).lower():
                    gx1 = idf_v600e
                    break
        if len(cholangiocarcinoma) != 0 and cholangiocarcinoma[0]:
            for eleK in KeywordsList:
                if 'braf' in str(eleK).lower():
                    gx = idf_braf
                    break
     # ---------------共现分析医学主题词-------------------#
        if len(cholangiocarcinoma) != 0 and cholangiocarcinoma[0]:
            for eleM in MeshHeadingNameList:
                if 'V600E' in eleM['MeshHeadingName']:
                    gx1 = idf_v600e
                    break
        if len(cholangiocarcinoma) != 0 and cholangiocarcinoma[0]:
            for eleM in MeshHeadingNameList:
                if 'B-raf' in eleM['MeshHeadingName']:
                    gx = idf_braf
                    break

        for ele in ChemicalNameList:
            if 'Cholangiocarcinoma' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_1
                break
        for ele in ChemicalNameList:
            if 'Bile Duct Neoplasms' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_2
                break
        for ele in ChemicalNameList:
            if 'Bile Ducts, Intrahepatic' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_3
                break
        for ele in ChemicalNameList:
            if 'Proto-Oncogene Proteins B-raf' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_4
                break
        for ele in ChemicalNameList:
            if 'Melanoma' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_5
                break
        for ele in ChemicalNameList:
            if 'Skin Neoplasms' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_6
                break


        for eleM in MeshHeadingNameList:
            if 'Cholangiocarcinoma' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_1
                break
        for eleM in MeshHeadingNameList:
            if 'Bile Duct Neoplasms' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_2
                break
        for eleM in MeshHeadingNameList:
            if 'Bile Ducts, Intrahepatic' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_3
                break
        for eleM in MeshHeadingNameList:
            if 'Proto-Oncogene Proteins B-raf' == eleM['MeshHeadingName']:
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
            if 'Aged' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_7
                break
        for eleM in MeshHeadingNameList:
            if 'Melanoma' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_8
                break
        for eleM in MeshHeadingNameList:
            if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_9
                break
        for eleM in MeshHeadingNameList:
            if 'Skin Neoplasms' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_10
                break


        for eleK in KeywordsList:
            if 'cholangiocarcinoma' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_1
                break
        for eleK in KeywordsList:
            if 'v600e' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_2
                break

        total_gx=gx1+gx2+gx3+gx
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
    sortsecond(mywords,mydata,6.5)
    count(mydata,mycount,"2")



