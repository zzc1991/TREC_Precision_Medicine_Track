import pymongo
import re
from math import log


myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["pubmed"]
mywords = mydb["freqwords1"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2019"]#pubmed中的主题词相关文献列表


mytopicdb=myclient["cs2019_bert"]
mydata=mytopicdb["cs2019_score_9"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2019_score_9_related"]#聚类后对应与主题相关联的文献


def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1 = 1.2
    b1 = 0.75
    k2 = 1.2
    b2 = 0.75
    idf_gastrointestinal = log((29138919 - 203456 + 0.5) / (203456 + 0.5), 10)
    idf_stromal = log((29138919 - 75278 + 0.5) / (75278 + 0.5), 10)
    idf_exon = log((29138919 - 62186 + 0.5) / (62186 + 0.5), 10)
    idf_kit = log((29138919 - 31493 + 0.5) / (31493 + 0.5), 10)
    idf_a502=log((29138919 - 61 + 0.5) / (61 + 0.5), 10)
    idf_502503 = log((29138919 - 12 + 0.5) / (12 + 0.5), 10)
    idf_ala502tyr503 = log((29138919 - 1 + 0.5) / (1 + 0.5), 10)

    idf_ele_1 = log((13670358 - 0 + 0.5) / (367 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 0 + 0.5) / (367 + 0.5), 10)
    idf_ele_4 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_5 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 5700 + 0.5) / (5700 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 46355 + 0.5) / (46355 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 5457 + 0.5) / (5457 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 11770 + 0.5) / (11770 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 8002162 + 0.5) / (8002162 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 2842020 + 0.5) / (2842020 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4029038 + 0.5) / (4029038 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 953 + 0.5) / (953 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 0 + 0.5) / (0 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        len_freq = 0
        gastrointestinal_score = 0
        stromal_score = 0
        exon_score = 0
        kit_score = 0
        a502_score=0
        s502503_score=0
        ala502tyr503_score=0
        gx = 0
        gx1 = 0
        gx2 = 0
        gx3 = 0
        cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
        ChemicalNameList = x['ChemicalNameList']
        MeshHeadingNameList = x['MeshHeadingNameList']
        KeywordsList = x['KeywordsList']
        wordfreq = x['wordfreq']

        gastrointestinal = [True for x in wordfreq.items() if 'gastrointestinal' in x]
        stromal = [True for x in wordfreq.items() if 'stromal' in x]
        kit = [True for x in wordfreq.items() if 'kit' in x]
        exon = [True for x in wordfreq.items() if 'exon' in x]
        a502 = [True for x in wordfreq.items() if 'a502' in x]
        y503 = [True for x in wordfreq.items() if 'y503' in x]
        s502503=[True for x in wordfreq.items() if '502-503' in x]
        ala502tyr503=[True for x in wordfreq.items() if 'ala502tyr503' in x]
        # ---------------摘要统计-------------------#
        for key in wordfreq:
            len_freq = len_freq + wordfreq[key]
        for key in wordfreq:
            if 'gastrointestinal' in key:
                gastrointestinal_score = gastrointestinal_score + wordfreq[key]
        for key in wordfreq:
            if 'stromal' in key:
                stromal_score = stromal_score + wordfreq[key]
        for key in wordfreq:
            if 'kit' in key:
                kit_score = kit_score + wordfreq[key]
        for key in wordfreq:
            if 'exon' in key:
                exon_score = exon_score + wordfreq[key]
        if len(a502) != 0 and a502[0] and len(y503) != 0 and y503[0]:
            for key in wordfreq:
                if 'a502' in key:
                    a502_score = a502_score + wordfreq[key]
        for key in wordfreq:
            if '502-503' in key:
                s502503_score = s502503_score + wordfreq[key]
        for key in wordfreq:
            if 'ala502tyr503' in key:
                ala502tyr503_score = ala502tyr503_score + wordfreq[key]

        #---------------共现分析摘要-------------------#
        if len(gastrointestinal) != 0 and gastrointestinal[0] and len(stromal) != 0 and stromal[0]:
            if len(kit) != 0 and kit[0] and len(exon) != 0 and exon[0]:
                gx=idf_kit+idf_exon
            if len(a502) != 0 and a502[0] and len(y503) != 0 and y503[0]:
                gx1=idf_a502
            if len(s502503) != 0 and s502503[0]:
                gx1=idf_502503
            if len(ala502tyr503) != 0 and ala502tyr503[0]:
                gx1=idf_ala502tyr503
        # ---------------共现分析化学-------------------#

        # ---------------共现分析医学主题词-------------------#

        # ---------------共现分析关键字-------------------#

        bm25_gastrointestinal_score = (((k1 + 1) * gastrointestinal_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + gastrointestinal_score))
        bm25_stromal_score = (((k1 + 1) * stromal_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + stromal_score))
        bm25_exon_score = (((k1 + 1) * exon_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + exon_score))
        bm25_kit_score = (((k1 + 1) * kit_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + kit_score))
        bm25_s502503_score = (((k1 + 1) * s502503_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + s502503_score))
        bm25_a502_score = (((k1 + 1) * a502_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + a502_score))
        bm25_ala502tyr503_score = (((k1 + 1) * ala502tyr503_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + ala502tyr503_score))

        bm25_ab_score = idf_gastrointestinal * bm25_gastrointestinal_score + idf_stromal * bm25_stromal_score
        +idf_exon * bm25_exon_score + idf_502503 * bm25_s502503_score+idf_kit*bm25_kit_score+idf_a502*bm25_a502_score+idf_ala502tyr503*bm25_ala502tyr503_score

        idf_para = [{str(gastrointestinal_score): idf_gastrointestinal}, {str(stromal_score): idf_stromal},
                    {str(exon_score): idf_exon}, {str(s502503_score): idf_502503},{str(kit_score): idf_kit},{str(ala502tyr503_score): idf_ala502tyr503},{str(a502_score): idf_a502}]

        for ele in ChemicalNameList:
            if 'Chromosome 9, trisomy' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_1
                break
        for ele in ChemicalNameList:
            if 'Trisomy' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_2
                break
        for ele in ChemicalNameList:
            if 'Chromosomes, Human, Pair 9' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_3
                break
        for ele in ChemicalNameList:
            if 'Exons' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_4
                break
        for ele in ChemicalNameList:
            if 'Gastrointestinal Stromal Tumors' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_5
                break

        for eleM in MeshHeadingNameList:
            if 'Gastrointestinal Stromal Tumors' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_1
                break
        for eleM in MeshHeadingNameList:
            if 'Exons' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_2
                break
        for eleM in MeshHeadingNameList:
            if 'Chromosomes, Human, Pair 9' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_3
                break
        for eleM in MeshHeadingNameList:
            if 'Trisomy' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_4
                break
        for eleM in MeshHeadingNameList:
            if 'Chromosome 9, trisomy' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_5
                break
        for eleM in MeshHeadingNameList:
            if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_6
                break
        for eleM in MeshHeadingNameList:
            if 'Male' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_7
                break
        for eleM in MeshHeadingNameList:
            if 'Middle Aged' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_9
                break
        for eleM in MeshHeadingNameList:
            if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_10
                break
        for eleK in KeywordsList:
            if 'gastrointestinal stromal tumor' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_1
                break
        for eleK in KeywordsList:
            if 'kit exon 9' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_2
                break

        total_gx = (gx + gx1 + gx2 + gx3)
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
    sortsecond(mywords,mydata,8)
    count(mydata,mycount,"9")


