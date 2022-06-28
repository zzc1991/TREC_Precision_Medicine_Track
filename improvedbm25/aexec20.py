from math import sqrt
import pymongo
import re
from math import  log

myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["bert"]
mywords = mydb["freqwords1"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2019"]#pubmed中的主题词相关文献列表


mytopicdb=myclient["cs2019_bert"]
mydata=mytopicdb["cs2019_score_20"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2019_score_20_related"]#聚类后对应与主题相关联的文献






def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1 = 1.2
    b1 = 0.75
    k2 = 1.2
    b2 = 0.75
    idf_pancreatic = log((29138919 - 167177 + 0.5) / (167177 + 0.5), 10)
    idf_brca2 = log((29138919 - 5591 + 0.5) / (5591 + 0.5), 10)

    idf_ele_1 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 1472 + 0.5) / (1472 + 0.5), 10)
    idf_ele_4 = log((13670358 - 3608+ 0.5) / (3608 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 67099 + 0.5) / (67099 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 12274 + 0.5) / (12274 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 3608 + 0.5) / (3608 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 8002162 + 0.5) / (8002162 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 4029038 + 0.5) / (4029038 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 260276 + 0.5) / (260276 + 0.5), 10)
    idf_eleM_11 = log((25389659 - 4029038 + 0.5) / (4029038 + 0.5), 10)
    idf_eleM_12 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 14780 + 0.5) / (14780 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 711 + 0.5) / (711 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        len_freq=0
        pancreatic_score = 0
        brca2_score = 0
        gx=0
        cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
        ChemicalNameList = x['ChemicalNameList']
        MeshHeadingNameList = x['MeshHeadingNameList']
        KeywordsList = x['KeywordsList']
        wordfreq = x['wordfreq']
        pancreatic = [True for x in wordfreq.items() if 'pancreatic' in x]
        # ---------------摘要统计-------------------#
        for key in wordfreq:
            len_freq = len_freq + wordfreq[key]


        for key in wordfreq:
            if 'pancreatic' in key:
                # ss3 = ss3 + 3
                pancreatic_score = pancreatic_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'brca2' == key1:
                # ss3 = ss3 + 3
                brca2_score = brca2_score + wordfreq[key]
        # ---------------共现分析摘要-------------------#
        if len(pancreatic) != 0 and pancreatic[0]:
                for key in wordfreq:
                    key = cop.sub('', key)
                    if 'brca2' == key:
                        gx=idf_brca2
                        break
    # ---------------共现分析化学-------------------#
        if len(pancreatic) != 0 and pancreatic[0]:
                for ele in ChemicalNameList:
                    if 'BRCA2' in ele['NameOfSubstance']:
                        gx = idf_brca2
                        break
    # ---------------共现分析关键字-------------------#
        if len(pancreatic) != 0 and pancreatic[0]:
                for eleK in KeywordsList:
                    if 'brca2' in str(eleK).lower():
                        gx = idf_brca2
                        break
        # ---------------共现分析医学主题词-------------------#
        if len(pancreatic) != 0 and pancreatic[0]:
                for eleM in MeshHeadingNameList:
                    if 'BRCA2' in eleM['MeshHeadingName']:
                        gx = idf_brca2
                        break

        bm25_pancreatic_score = (((k1 + 1) * pancreatic_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + pancreatic_score))
        bm25_brca2_score = (((k1 + 1) * brca2_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + brca2_score))

        bm25_ab_score = idf_pancreatic * bm25_pancreatic_score + idf_brca2 * bm25_brca2_score

        idf_para = [{str(pancreatic_score): idf_pancreatic}, {str(brca2_score): idf_brca2}]
        for ele in ChemicalNameList:
            # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
            if 'Pancreatic Neoplasms' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_1
                break
        for ele in ChemicalNameList:
            # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
            if 'Pancreatectomy' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_2
                break
        for ele in ChemicalNameList:
            # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
            if 'BRCA2 protein, human' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_3
                break
        for ele in ChemicalNameList:
            # if re.findall(r'(BRAF|Proto-Oncogene Proteins B-raf|human|humans|male)',ele['NameOfSubstance']):
            if 'BRCA2 Protein' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_4
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'Pancreatic Neoplasms' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_1
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'Pancreatectomy' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_2
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'BRCA2 protein, human' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_3
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'BRCA2 Protein' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_4
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_6
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'Male' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_7
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'Middle Aged' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_8
                break
        for eleM in MeshHeadingNameList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_9
                break
        for eleK in KeywordsList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'pancreatic' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_1
                break
        for eleK in KeywordsList:
            # if re.findall(r'(Melanoma|Proto-Oncogene Proteins B-raf|Humans|Neoplasms|Neoplasm|Male|Mutation|Mutational)',eleM['MeshHeadingName']):
            if 'brca2' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_2
                break

        total_gx = gx
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
    sortsecond(mywords,mydata,6)
    count(mydata,mycount,"20")



