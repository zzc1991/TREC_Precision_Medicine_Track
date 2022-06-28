import pymongo
import re
from math import log


myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["bert"]
mywords = mydb["freqwords1"] #pubmed中所有的词频、化学词、关键词和主题词表
mytopic=mydb["topics2019"]#pubmed中的主题词相关文献列表


mytopicdb=myclient["cs2019_bert"]
mydata=mytopicdb["cs2019_score_14"]#按词表长度改进过后的2次排序表
mycount = mytopicdb["cs2019_score_14_related"]#聚类后对应与主题相关联的文献





def sortsecond(myfreq,mydata,yuzhi):
    k = 0
    k1 = 1.2
    b1 = 0.75
    k2 = 1.2
    b2 = 0.75
    idf_colon = log((29138919 - 133274 + 0.5) / (133274 + 0.5), 10)
    idf_mlh1 = log((29138919 - 4948 + 0.5) / (4948 + 0.5), 10)
    idf_methylation = log((29138919 - 99607 + 0.5) / (99607 + 0.5), 10)
    idf_suppression = log((29138919 - 276766 + 0.5) / (276766 + 0.5), 10)
    idf_microsatellite = log((29138919 - 38274 + 0.5) / (38274 + 0.5), 10)
    idf_instability = log((29138919 - 100164 + 0.5) / (100164 + 0.5), 10)

    idf_ele_1 = log((13670358 - 0 + 0.5) / (367 + 0.5), 10)
    idf_ele_2 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_3 = log((13670358 - 3004 + 0.5) / (3004 + 0.5), 10)
    idf_ele_4 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)
    idf_ele_5 = log((13670358 - 0 + 0.5) / (0 + 0.5), 10)

    idf_eleM_1 = log((25389659 - 67366 + 0.5) / (67366 + 0.5), 10)
    idf_eleM_2 = log((25389659 - 2662+ 0.5) / (2662 + 0.5), 10)
    idf_eleM_3 = log((25389659 - 3004 + 0.5) / (3004 + 0.5), 10)
    idf_eleM_4 = log((25389659 - 42230 + 0.5) / (42230 + 0.5), 10)
    idf_eleM_5 = log((25389659 - 17437618 + 0.5) / (17437618 + 0.5), 10)
    idf_eleM_6 = log((25389659 - 8002162 + 0.5) / (8002162 + 0.5), 10)
    idf_eleM_7 = log((25389659 - 2842020 + 0.5) / (2842020 + 0.5), 10)
    idf_eleM_8 = log((25389659 - 0 + 0.5) / (0 + 0.5), 10)
    idf_eleM_9 = log((25389659 - 4785026 + 0.5) / (4785026 + 0.5), 10)

    idf_eleK_1 = log((5435471 - 22646 + 0.5) / (22646 + 0.5), 10)
    idf_eleK_2 = log((5435471 - 200 + 0.5) / (200 + 0.5), 10)
    idf_eleK_3 = log((5435471 - 7 + 0.5) / (7 + 0.5), 10)
    for x in myfreq.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):
        ss1 = 0
        ss2 = 0
        ss4 = 0
        len_freq = 0
        colon_score = 0
        mlh1_score = 0
        methylation_score = 0
        suppression_score = 0
        microsatellite_score = 0
        instability_score = 0
        gx = 0
        gx1 = 0
        gx2 = 0
        gx3 = 0
        gx4=0
        cop = re.compile("[^\u4e00-\u9fa5^a-z^A-Z^0-9]")  # 匹配不是中文、大小写、数字的其他字符
        ChemicalNameList = x['ChemicalNameList']
        MeshHeadingNameList = x['MeshHeadingNameList']
        KeywordsList = x['KeywordsList']
        wordfreq = x['wordfreq']

        colon = [True for x in wordfreq.items() if 'colon' in x]
        # ---------------摘要统计-------------------#
        for key in wordfreq:
            len_freq = len_freq + wordfreq[key]

        for key in wordfreq:
            if 'colon' in key:
                colon_score = colon_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'mlh1' == key1:
                mlh1_score = mlh1_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'methylation' in key1:
                methylation_score = methylation_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'suppression' in key1:
                suppression_score = suppression_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'microsatellite' in key1:
                microsatellite_score = microsatellite_score + wordfreq[key]
        for key in wordfreq:
            key1 = cop.sub('', key)
            if 'instability' in key1:
                instability_score = instability_score + wordfreq[key]

        #---------------共现分析摘要-------------------#
        if len(colon) != 0 and colon[0]:
            for key in wordfreq:
                if 'mlh1' == key:
                    gx = idf_mlh1
                if 'methylation' in key:
                    gx1 = idf_methylation
                if 'suppression' in key:
                    gx2 = idf_suppression
                if 'microsatellite' in key:
                    gx3 = idf_microsatellite
                if 'instability' in key:
                    gx4 = idf_instability
                    break
        # ---------------共现分析化学-------------------#

        if len(colon) != 0 and colon[0]:
            for ele in ChemicalNameList:
                if 'MLH1' in ele['NameOfSubstance']:
                    gx = idf_mlh1
                    break

        # ---------------共现分析医学主题词-------------------#


        if len(colon) != 0 and colon[0]:
            for eleM in MeshHeadingNameList:
                if 'MLH1' in eleM['MeshHeadingName']:
                    gx = idf_mlh1
                    break



        # ---------------共现分析关键字-------------------#

        if len(colon) != 0 and colon[0]:
            for eleK in KeywordsList:
                if 'mlh1' in str(eleK).lower():
                    gx = idf_mlh1
                    break

        bm25_colon_score = (((k1 + 1) * colon_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + colon_score))
        bm25_mlh1_score = (((k1 + 1) * mlh1_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + mlh1_score))
        bm25_methylation_score = (((k1 + 1) * methylation_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + methylation_score))
        bm25_suppression_score = (((k1 + 1) * suppression_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + suppression_score))
        bm25_microsatellite_score = (((k1 + 1) * microsatellite_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + microsatellite_score))
        bm25_instability_score = (((k1 + 1) * instability_score) / ((k1 * (b1 + (1 - b1) * (len_freq / 83))) + instability_score))

        bm25_ab_score = idf_colon * bm25_colon_score + idf_mlh1 * bm25_mlh1_score + idf_methylation * bm25_methylation_score + idf_suppression * bm25_suppression_score
        +idf_microsatellite*bm25_microsatellite_score+idf_instability*bm25_instability_score

        idf_para = [{str(colon_score): idf_colon}, {str(mlh1_score): idf_mlh1},
                    {str(methylation_score): idf_methylation}, {str(suppression_score): idf_suppression}, {str(microsatellite_score): idf_microsatellite}, {str(instability_score): idf_instability}]

        for ele in ChemicalNameList:
            if 'Colonic Neoplasms' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_1
                break
        for ele in ChemicalNameList:
            if 'Microsatellite Instability' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_2
                break
        for ele in ChemicalNameList:
            if 'MutL Protein Homolog 1' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_3
                break
        for ele in ChemicalNameList:
            if 'DNA Methylation' == ele['NameOfSubstance']:
                ss1 = ss1 + idf_ele_4
                break


        for eleM in MeshHeadingNameList:
            if 'Colonic Neoplasms' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_1
                break
        for eleM in MeshHeadingNameList:
            if 'Microsatellite Instability' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_2
                break
        for eleM in MeshHeadingNameList:
            if 'MutL Protein Homolog 1' == eleM['MeshHeadingName']:
                ss2 = ss2 + idf_eleM_3
                break

        for eleM in MeshHeadingNameList:
            if 'DNA Methylation' == eleM['MeshHeadingName']:
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
            if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                ss2 = ss2 + idf_eleM_9
                break

        for eleK in KeywordsList:
            if 'colon' in str(eleK).lower():
                ss4 = ss4 +idf_eleK_1
                break
        for eleK in KeywordsList:
            if 'mlh1' in str(eleK).lower():
                ss4 = ss4 + idf_eleK_2
                break
        total_gx = (gx + gx1 + gx2 + gx3+gx4)
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
    count(mydata,mycount,"14")


