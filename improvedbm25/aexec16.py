import pymongo
import re
from math import log


myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["bert"]
mywords = mydb["freqwords1"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2019"]#pubmed中的主题词相关文献列表


mytopicdb=myclient["cs2019_bert"]
mydata=mytopicdb["cs2019_score_16"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2019_score_16_related"]#聚类后对应与主题相关联的文献


def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1=1.2
    k2=1.2
    b1=0.75
    b2=0.75
    idf_papillary = log((29138919 - 50419 + 0.5) / (50419 + 0.5), 10)
    idf_thyroid = log((29138919 - 161588 + 0.5) / (161588 + 0.5), 10)
    idf_carcinoma = log((29138919 - 494907 + 0.5) / (494907 + 0.5), 10)
    idf_braf = log((29138919 - 8527 + 0.5) / (8527 + 0.5), 10)
    idf_v600e = log((29138919 - 1887 + 0.5) / (1887 + 0.5), 10)


    idf_ele_1 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_4 = log((13670358 - 7222 + 0.5) / (7222 + 0.5), 10)
    idf_ele_5 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_6 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 3212 + 0.5) / (3212 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 46290 + 0.5) / (46290 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 18540 + 0.5) / (18540 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 7221 + 0.5) / (7221 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 8104149 + 0.5) / (8104149 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 4029038 + 0.5) / (4029038 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 46290 + 0.5) / (46290 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 396072 + 0.5) / (396072 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 900 + 0.5) / (900 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 415 + 0.5) / (415 + 0.5), 10)
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
        papillary_score=0
        thyroid_score = 0
        carcinoma_score=0
        braf_score = 0
        v600e_score = 0
        if int(x['PMID']) <= 27868941:
                cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
                ChemicalNameList = x['ChemicalNameList']
                MeshHeadingNameList = x['MeshHeadingNameList']
                KeywordsList = x['KeywordsList']
                wordfreq = x['wordfreq']
                papillary = [True for x in wordfreq.items() if 'papillary' in x]
                thyroid = [True for x in wordfreq.items() if 'thyroid' in x]
                carcinoma = [True for x in wordfreq.items() if 'carcinoma' in x]


                # ---------------摘要统计-------------------#


                for key in wordfreq:
                    len_freq = len_freq + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'papillary' in key1:
                        papillary_score = papillary_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'thyroid' in key1:
                        thyroid_score = thyroid_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'carcinoma' in key1:
                        carcinoma_score = carcinoma_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'braf' in key1:
                        braf_score = braf_score + wordfreq[key]
                for key in wordfreq:
                    key1 = cop.sub('', key)
                    if 'v600e' in key1:
                        v600e_score = v600e_score + wordfreq[key]


                bm25_papillary_score = (((k1+1)*papillary_score)/((k1*(b1+(1-b1)*(len_freq/85)))+papillary_score))
                bm25_thyroid_score = (((k1 + 1) * thyroid_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + thyroid_score))
                bm25_carcinoma_score = (((k1 + 1) * carcinoma_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + carcinoma_score))
                bm25_braf_score = (((k1 + 1) * braf_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + braf_score))
                bm25_v600e_score = (((k1 + 1) * v600e_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 85))) + v600e_score))


                bm25_ab_score =idf_papillary*bm25_papillary_score+idf_thyroid*bm25_thyroid_score+idf_carcinoma*bm25_carcinoma_score+idf_braf*bm25_braf_score+idf_v600e*bm25_v600e_score

                idf_para=[{str(papillary_score):idf_papillary},{str(thyroid_score):idf_thyroid},{str(carcinoma_score):idf_carcinoma},{str(braf_score):idf_braf},{str(v600e_score):idf_v600e}]

                # ---------------共现分析摘要-------------------#
                if len(papillary) != 0 and papillary[0] and len(thyroid) != 0 and thyroid[0] and len(carcinoma) != 0 and carcinoma[0]:
                    for key in wordfreq:
                        key = cop.sub('', key)
                        if 'braf' in key:
                            gx = idf_braf
                        if 'v600e' in key:
                            gx1 = idf_v600e

            # ---------------共现分析化学-------------------#
                if len(papillary) != 0 and papillary[0] and len(thyroid) != 0 and thyroid[0] and len(carcinoma) != 0 and carcinoma[0]:
                    for ele in ChemicalNameList:
                        if 'V600E' in ele['NameOfSubstance']:
                            gx1 = idf_v600e
                            break
            # ---------------共现分析关键字-------------------#
                if len(papillary) != 0 and papillary[0] and len(thyroid) != 0 and thyroid[0] and len(carcinoma) != 0 and carcinoma[0]:
                    for eleK in KeywordsList:
                        if 'v600e' in str(eleK).lower():
                            gx1= idf_v600e
                            break
             # ---------------共现分析医学主题词-------------------#
                if len(papillary) != 0 and papillary[0] and len(thyroid) != 0 and thyroid[0] and len(carcinoma) != 0 and carcinoma[0]:
                    for eleM in MeshHeadingNameList:
                        if 'V600E' in eleM['MeshHeadingName']:
                            gx1 = idf_v600e
                            break

                for ele in ChemicalNameList:
                    if 'Thyroid Cancer, Papillary' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_1
                        break
                for ele in ChemicalNameList:
                    if 'Thyroid Neoplasms' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_2
                        break
                for ele in ChemicalNameList:
                    if 'Carcinoma, Papillary' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_3
                        break
                for ele in ChemicalNameList:
                    if 'Proto-Oncogene Proteins B-raf' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_4
                        break
                for ele in ChemicalNameList:
                    if 'Thyroid Neoplasms' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_5
                        break
                for ele in ChemicalNameList:
                    if 'Mutation' in ele['NameOfSubstance']:
                        ss1 = ss1 + idf_ele_6
                        break

                for eleM in MeshHeadingNameList:
                    if 'Thyroid Cancer, Papillary' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_1
                        break
                for eleM in MeshHeadingNameList:
                    if 'Thyroid Neoplasms' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_2
                        break
                for eleM in MeshHeadingNameList:
                    if 'Carcinoma, Papillary' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_3
                        break
                for eleM in MeshHeadingNameList:
                    if 'Proto-Oncogene Proteins B-raf' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_4
                        break
                for eleM in MeshHeadingNameList:
                    if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_5
                        break
                for eleM in MeshHeadingNameList:
                    if 'Female' in eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_6
                        break
                for eleM in MeshHeadingNameList:
                    if 'Middle Aged' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_7
                        break
                for eleM in MeshHeadingNameList:
                    if 'Thyroid Neoplasms' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_8
                        break
                for eleM in MeshHeadingNameList:
                    if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                        ss2 = ss2 + idf_eleM_9
                        break
                for eleM in MeshHeadingNameList:
                    if 'Mutation' == eleM['MeshHeadingName']:
                        ss2 = ss2 + idf_eleM_10
                        break

                for eleK in KeywordsList:
                    if 'papillary thyroid carcinoma' in str(eleK).lower():
                        ss4 = ss4 + idf_eleK_1
                        break
                for eleK in KeywordsList:
                    if 'v600e' in str(eleK).lower():
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
    sortsecond(mywords,mydata,13)
    count(mydata,mycount,"16")



