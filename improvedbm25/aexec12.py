import pymongo
import re
from math import log


myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["bert"]
mywords = mydb["freqwords1"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2019"]#pubmed中的主题词相关文献列表


mytopicdb=myclient["cs2019_bert"]
mydata=mytopicdb["cs2019_score_12"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2019_score_12_related"]#聚类后对应与主题相关联的文献


def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1 = 1.2
    b1 = 0.75
    k2 = 1.2
    b2 = 0.75
    idf_inflammatory = log((29138919 - 645801 + 0.5) / (645801 + 0.5), 10)
    idf_myofibroblastic = log((29138919 - 3440 + 0.5) / (3440 + 0.5), 10)
    idf_ranbp2alk = log((29138919 - 11 + 0.5) / (11 + 0.5), 10)
    idf_fusion = log((29138919 - 167225 + 0.5) / (167225 + 0.5), 10)

    idf_ele_1 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 212 + 0.5) / (212 + 0.5), 10)
    idf_ele_4 = log((13670358 - 2505 + 0.5) / (2505 + 0.5), 10)
    idf_ele_5 = log((13670358 - 14529 + 0.5) / (14529 + 0.5), 10)
    idf_ele_6 = log((13670358 - 16271 + 0.5) / (16271 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 2303 + 0.5) / (2303 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 391874 + 0.5) / (391874 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 2505 + 0.5) / (2505 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 14529 + 0.5) / (14529 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 8104149 + 0.5) / (8104149 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 16271 + 0.5) / (16271 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 702001 + 0.5) / (702001 + 0.5), 10)
    idf_eleM_10 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 225 + 0.5) / (225 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 5 + 0.5) / (5 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        len_freq = 0
        inflammatory_score = 0
        myofibroblastic_score = 0
        ranbp2alk_score = 0
        fusion_score = 0

        gx = 0
        gx1 = 0
        gx2 = 0
        gx3 = 0
        cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
        ChemicalNameList = x['ChemicalNameList']
        MeshHeadingNameList = x['MeshHeadingNameList']
        KeywordsList = x['KeywordsList']
        wordfreq = x['wordfreq']

        inflammatory = [True for x in wordfreq.items() if 'inflammatory' in x]
        myofibroblastic = [True for x in wordfreq.items() if 'myofibroblastic' in x]
        # ---------------摘要统计-------------------#
        for key in wordfreq:
            len_freq = len_freq + wordfreq[key]
        for key in wordfreq:
            if 'inflammatory' in key:
                inflammatory_score = inflammatory_score + wordfreq[key]
        for key in wordfreq:
            if 'myofibroblastic' in key:
                myofibroblastic_score = myofibroblastic_score + wordfreq[key]
        for key in wordfreq:
            if 'ranbp2-alk' in key:
                ranbp2alk_score = ranbp2alk_score + wordfreq[key]
        for key in wordfreq:
            if 'fusion' in key:
                fusion_score = fusion_score + wordfreq[key]

        #---------------共现分析摘要-------------------#
        if len(inflammatory) != 0 and inflammatory[0] and len(myofibroblastic) != 0 and myofibroblastic[0]:
            for key in wordfreq:
                if 'ranbp2-alk' in key:
                    gx = idf_ranbp2alk
                if 'fusion' in key:
                    gx1 = idf_fusion
        # ---------------共现分析化学-------------------#
        if len(inflammatory) != 0 and inflammatory[0] and len(myofibroblastic) != 0 and myofibroblastic[0]:
            for ele in ChemicalNameList:
                if 'RANBP2-ALK' in ele['NameOfSubstance']:
                    gx = idf_ranbp2alk

        # ---------------共现分析医学主题词-------------------#
        if len(inflammatory) != 0 and inflammatory[0] and len(myofibroblastic) != 0 and myofibroblastic[0]:
            for eleM in MeshHeadingNameList:
                if 'RANBP2-ALK' in eleM['MeshHeadingName']:
                    gx = idf_ranbp2alk
        # ---------------共现分析关键字-------------------#
        if len(inflammatory) != 0 and inflammatory[0] and len(myofibroblastic) != 0 and myofibroblastic[0]:
            for eleK in KeywordsList:
                if 'ranbp2-alk' in str(eleK).lower():
                    gx = idf_ranbp2alk

        bm25_inflammatory_score = (((k1 + 1) * inflammatory_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + inflammatory_score))
        bm25_myofibroblastic_score = (((k1 + 1) * myofibroblastic_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + myofibroblastic_score))
        bm25_ranbp2alk_score = (((k1 + 1) * ranbp2alk_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + ranbp2alk_score))
        bm25_fusion_score = (((k1 + 1) * fusion_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + fusion_score))

        bm25_ab_score = idf_inflammatory * bm25_inflammatory_score + idf_myofibroblastic * bm25_myofibroblastic_score+idf_ranbp2alk*bm25_ranbp2alk_score+idf_fusion*bm25_fusion_score

        idf_para = [{str(inflammatory_score): idf_inflammatory}, {str(myofibroblastic_score): idf_myofibroblastic},
                    {str(ranbp2alk_score): idf_ranbp2alk}, {str(fusion_score): idf_fusion}]

        for ele in ChemicalNameList:
            if 'Granuloma, Plasma Cell' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_1
                break
        for ele in ChemicalNameList:
            if 'Neoplasms' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_2
                break
        for ele in ChemicalNameList:
            if 'ran-binding protein 2' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_3
                break
        for ele in ChemicalNameList:
            if 'Nuclear Pore Complex Proteins' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_4
                break
        for ele in ChemicalNameList:
            if 'Molecular Chaperones' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_5
                break
        for ele in ChemicalNameList:
            if 'Receptor Protein-Tyrosine Kinases' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_6
                break

        for eleM in MeshHeadingNameList:
            if 'Granuloma, Plasma Cell' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_1
                break
        for eleM in MeshHeadingNameList:
            if 'Neoplasms' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_2
                break
        for eleM in MeshHeadingNameList:
            if 'ran-binding protein 2' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_3
                break
        for eleM in MeshHeadingNameList:
            if 'Nuclear Pore Complex Proteins' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_4
                break
        for eleM in MeshHeadingNameList:
            if 'Molecular Chaperones' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_5
                break
        for eleM in MeshHeadingNameList:
            if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_6
                break
        for eleM in MeshHeadingNameList:
            if 'Female' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_7
                break
        for eleM in MeshHeadingNameList:
            if 'Receptor Protein-Tyrosine Kinases' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_8
                break
        for eleM in MeshHeadingNameList:
            if 'Young' in eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_9
                break
        for eleM in MeshHeadingNameList:
            if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_10
                break
        for eleK in KeywordsList:
            if 'inflammatory myofibroblastic tumor' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_1
                break
        for eleK in KeywordsList:
            if 'ranbp2-alk' in str(eleK).lower():
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
    count(mydata,mycount,"12")


