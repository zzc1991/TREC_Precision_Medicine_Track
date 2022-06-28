import pymongo
import re
from math import log


myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["bert"]
mywords = mydb["freqwords1"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2019"]#pubmed中的主题词相关文献列表


mytopicdb=myclient["cs2019_bert"]
mydata=mytopicdb["cs2019_score_21"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2019_score_21_related"]#聚类后对应与主题相关联的文献


def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1 = 1.2
    b1 = 0.75
    k2 = 1.2
    b2 = 0.75
    idf_gastrointestinal = log((29138919 - 203456 + 0.5) / (203456 + 0.5), 10)
    idf_stromal = log((29138919 - 75278 + 0.5) / (75278 + 0.5), 10)
    idf_kit = log((29138919 - 31493 + 0.5) / (31493 + 0.5), 10)
    idf_v654a = log((29138919 - 18 + 0.5) / (18 + 0.5), 10)

    idf_ele_1 = log((13670358 - 0 + 0.5) / (367 + 0.5), 10)
    idf_ele_2 = log((13670358 - 9503 + 0.5) / (9503 + 0.5), 10)
    idf_ele_3 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_4 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_5 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 5700 + 0.5) / (5700 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 9503 + 0.5) / (9503 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 8002162 + 0.5) / (8002162 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 2842020 + 0.5) / (2842020 + 0.5), 10)
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
        kit_score = 0
        v654a_score = 0

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
            if 'v654a' in key:
                v654a_score = v654a_score + wordfreq[key]

        #---------------共现分析摘要-------------------#
        if len(gastrointestinal) != 0 and gastrointestinal[0] and len(stromal) != 0 and stromal[0]:
            for key in wordfreq:
                key = cop.sub('', key)
                if 'kit' in key:
                    gx = idf_kit
                if 'v654a' in key:
                    gx1 = idf_v654a
        # ---------------共现分析化学-------------------#
        if len(gastrointestinal) != 0 and gastrointestinal[0] and len(stromal) != 0 and stromal[0]:
            for ele in ChemicalNameList:
                if 'V654A' in ele['NameOfSubstance']:
                    gx1 = idf_v654a
                    break
        # ---------------共现分析医学主题词-------------------#
        if len(gastrointestinal) != 0 and gastrointestinal[0] and len(stromal) != 0 and stromal[0]:
            for eleM in MeshHeadingNameList:
                if 'V654A' in eleM['MeshHeadingName']:
                    gx1 = idf_v654a
                    break
        # ---------------共现分析关键字-------------------#
        if len(gastrointestinal) != 0 and gastrointestinal[0] and len(stromal) != 0 and stromal[0]:
            for eleK in KeywordsList:
                if 'v654a' in str(eleK).lower():
                    gx1 = idf_v654a
                    break

        bm25_gastrointestinal_score = (((k1 + 1) * gastrointestinal_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + gastrointestinal_score))
        bm25_stromal_score = (((k1 + 1) * stromal_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + stromal_score))
        bm25_kit_score = (((k1 + 1) * kit_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + kit_score))
        bm25_v654a_score = (((k1 + 1) * v654a_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + v654a_score))

        bm25_ab_score = idf_gastrointestinal * bm25_gastrointestinal_score + idf_stromal * bm25_stromal_score+idf_kit*bm25_kit_score+idf_v654a*bm25_v654a_score

        idf_para = [{str(gastrointestinal_score): idf_gastrointestinal}, {str(stromal_score): idf_stromal},
                    {str(kit_score): idf_kit}, {str(v654a_score): idf_v654a}]

        for ele in ChemicalNameList:
            if 'Gastrointestinal Stromal Tumors' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_1
                break
        for ele in ChemicalNameList:
            if 'Imatinib Mesylate' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_2
                break



        for eleM in MeshHeadingNameList:
            if 'Gastrointestinal Stromal Tumors' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_1
                break
        for eleM in MeshHeadingNameList:
            if 'Imatinib Mesylate' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_2
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
            if 'Aged' == eleM['MeshHeadingName']:
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
            if 'v654a' in str(eleK).lower():
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
    count(mydata,mycount,"21")


