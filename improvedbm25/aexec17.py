import pymongo
import re
from math import log

myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["bert"]
mywords = mydb["freqwords1"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2019"]#pubmed中的主题词相关文献列表


mytopicdb=myclient["cs2019_bert"]
mydata=mytopicdb["cs2019_score_17"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2019_score_17_related"]#聚类后对应与主题相关联的文献






def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1 = 1.2
    b1 = 0.75
    k2 = 1.2
    b2 = 0.75
    idf_glioblastoma = log((29138919 - 28553 + 0.5) / (28553 + 0.5), 10)
    idf_idh1 = log((29138919 - 1470 + 0.5) / (1470 + 0.5), 10)
    idf_r132h = log((29138919 - 217 + 0.5) / (217 + 0.5), 10)

    idf_ele_1 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_4 = log((13670358 - 5676 + 0.5) / (5676 + 0.5), 10)
    idf_ele_5 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 22078 + 0.5) / (22078 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 103237 + 0.5) / (103237 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 255959 + 0.5) / (255959 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 5676 + 0.5) / (5676 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 5676 + 0.5) / (5676 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 36507 + 0.5) / (36507 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 7716 + 0.5) / (7716 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 8002162 + 0.5) / (8002162 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 4029038 + 0.5) / (4029038 + 0.5), 10)
    idf_eleM_11 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)
    idf_eleM_12 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 5823 + 0.5) / (5823 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 26 + 0.5) / (26 + 0.5), 10)
    idf_eleK_3 = log((5435471 - 333 + 0.5) / (333 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        len_freq = 0
        glioblastoma_score = 0
        idh1_score = 0
        r132h_score = 0
        gx = 0
        gx1 = 0
        gx2 = 0
        gx3 = 0
        gx4=0
        gx5=0
        cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
        ChemicalNameList = x['ChemicalNameList']
        MeshHeadingNameList = x['MeshHeadingNameList']
        KeywordsList = x['KeywordsList']
        wordfreq = x['wordfreq']
        glioblastoma = [True for x in wordfreq.items() if 'glioblastoma' in x]
        # ---------------摘要统计-------------------#
        for key in wordfreq:
            len_freq = len_freq + wordfreq[key]
        for key in wordfreq:
            if 'glioblastoma' in key:
                # ss3 = ss3 + 3
                glioblastoma_score = glioblastoma_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'idh1' in key1:
                # ss3 = ss3 + 3
                idh1_score = idh1_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'r132h' in key1:
                # ss3 = ss3 + 3
                r132h_score = r132h_score + wordfreq[key]

        # ---------------共现分析摘要-------------------#
        if len(glioblastoma) != 0 and glioblastoma[0]:
                for key in wordfreq:
                    key = cop.sub('', key)
                    if 'r132h' in key:
                        gx=idf_r132h
                        break

        if len(glioblastoma) != 0 and glioblastoma[0]:
                for key in wordfreq:
                    key = cop.sub('', key)
                    if 'idh1' in key:
                        gx1=idf_idh1
                        break


    # ---------------共现分析化学-------------------#
        if len(glioblastoma) != 0 and glioblastoma[0]:
                for ele in ChemicalNameList:
                    if 'R132H' in ele['NameOfSubstance']:
                        gx = idf_r132h
                        break
        if len(glioblastoma) != 0 and glioblastoma[0]:
                for ele in ChemicalNameList:
                    if 'IDH1' in ele['NameOfSubstance']:
                        gx1 = idf_idh1
                        break

    # ---------------共现分析关键字-------------------#
        if len(glioblastoma) != 0 and glioblastoma[0]:
                for eleK in KeywordsList:
                    if 'r132h' in str(eleK).lower():
                        gx = idf_r132h
                        break

        if len(glioblastoma) != 0 and glioblastoma[0]:
                for eleK in KeywordsList:
                    if 'idh1' in str(eleK).lower():
                        gx1 = idf_idh1
                        break

        # ---------------共现分析医学主题词-------------------#
        for key in wordfreq:
            if 'glioblastoma' in key:
                for eleM in MeshHeadingNameList:
                    if 'R132H' in eleM['MeshHeadingName']:
                        gx = idf_r132h
                        break

        for key in wordfreq:
            if 'glioblastoma' in key:
                for eleM in MeshHeadingNameList:
                    if 'IDH1' in eleM['MeshHeadingName']:
                        gx1 = idf_idh1
                        break
        bm25_glioblastoma_score = (((k1 + 1) * glioblastoma_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + glioblastoma_score))
        bm25_idh1_score = (((k1 + 1) * idh1_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + idh1_score))
        bm25_r132h_score = (((k1 + 1) * r132h_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + r132h_score))

        bm25_ab_score = idf_glioblastoma * bm25_glioblastoma_score + idf_idh1 * bm25_idh1_score + idf_r132h * bm25_r132h_score

        idf_para = [{str(glioblastoma_score): idf_glioblastoma}, {str(idh1_score): idf_idh1},
                    {str(r132h_score): idf_r132h}]
        for ele in ChemicalNameList:
            # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
            if 'Glioblastoma' in ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_1
                break
        for ele in ChemicalNameList:
            # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
            if 'Brain Neoplasms' in ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_2
                break
        for ele in ChemicalNameList:
            # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
            if 'Cell Line, Tumor' in ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_3
                break
        for ele in ChemicalNameList:
            # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
            if 'Isocitrate Dehydrogenase' in ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_4
                break
        for ele in ChemicalNameList:
            # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
            if 'Glioma' in ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_5
                break

        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'Glioblastoma' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_1
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'Brain Neoplasms' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_2
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'Cell Line, Tumor' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_3
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'Isocitrate Dehydrogenase' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_5
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'Glioma' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_6
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_8
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'Male' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_9
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'Middle Aged' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_10
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_11
                break

        for eleK in KeywordsList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'glioblastoma' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_1
                break
        for eleK in KeywordsList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'r132h' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_2
                break
        for eleK in KeywordsList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'idh1' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_3
                break

        total_gx = (gx + gx1 + gx2 + gx3+gx4+gx5)
        cmk_len = len(ChemicalNameList) + len(MeshHeadingNameList) + len(KeywordsList)
        bm25_cmk_len = ss1 + ss2 + ss4
        bm25_cmk_score = (((k2 + 1) * bm25_cmk_len) / ((k2 * (b2 + (1 - b2) * (cmk_len / 13))) + bm25_cmk_len))
        bm25_score = bm25_ab_score + bm25_cmk_score + total_gx

        if (bm25_score > yuzhi):
            mydict = {"PMID": x['PMID'], "ab_score": bm25_ab_score, "idf_para": idf_para,
                      "cmk_len": cmk_len, "cmk_freq": bm25_cmk_len, "bm25_cmk_score": bm25_cmk_score,
                      "gx": total_gx,
                      "bm25_score": bm25_score,
                      "ChemicalNameList": x['ChemicalNameList'],
                      "MeshHeadingNameList": x['MeshHeadingNameList'], "KeywordsList": x['KeywordsList']}
            y = mydata.insert_one(mydict)
            k = k + 1
            print(str(y) + '---------' + str(k))

def count(mysort,mycount,topic):
    for x in mysort.find({},
                         {'PMID', 'ab_score', 'idf_para', 'cmk_len', 'cmk_freq', 'bm25_cmk_score', 'gx', 'bm25_score',
                          'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'}):
        kk = 0
        for y in mytopic.find({"topic": topic}, {'PMID', 'relate'}):
            if x['PMID'] == y['PMID']:
                mydict = {"PMID": x['PMID'], "related": y['relate'], "ab_score": x["ab_score"],
                          "idf_para": x['idf_para'],
                          "cmk_len": x['cmk_len'], "cmk_freq": x['cmk_freq'], 'bm25_cmk_score': x['bm25_cmk_score'],
                          'gx': x['gx'],
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
    count(mydata,mycount,"17")



