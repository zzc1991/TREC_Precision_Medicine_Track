import pymongo
import re

myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mydb = myclient["pubmed"]
mywords = mydb["freqwords3"] #pubmed中所有的词频、化学词、关键词和主题词表

if __name__ == '__main__':
    i=0
    idf_ele_1 = 0
    idf_ele_2 = 0
    idf_ele_3 = 0
    idf_ele_4 = 0
    idf_ele_5 = 0
    idf_ele_6 = 0
    idf_ele_7 = 0
    idf_ele_8 = 0

    idf_eleM_1 = 0
    idf_eleM_2 = 0
    idf_eleM_3 = 0
    idf_eleM_4 = 0
    idf_eleM_5 = 0
    idf_eleM_6 = 0
    idf_eleM_7 = 0
    idf_eleM_8 = 0
    idf_eleM_9 = 0
    idf_eleM_10 = 0
    idf_eleM_11 = 0
    idf_eleM_12 = 0

    idf_eleK_1 = 0
    idf_eleK_2 = 0
    idf_eleK_3 = 0

    squamous=0
    for x in mywords.find({}, {'PMID', 'wordfreq', 'ChemicalNameList', 'MeshHeadingNameList', 'KeywordsList'},
                         no_cursor_timeout=True):

        ChemicalNameList = x['ChemicalNameList']
        MeshHeadingNameList = x['MeshHeadingNameList']
        KeywordsList = x['KeywordsList']
        for ele in ChemicalNameList:
            if 'Malignant Hyperthermia' == ele['NameOfSubstance']:
                idf_ele_1 = idf_ele_1 + 1
                break
        for ele in ChemicalNameList:
            if 'Muscle Contraction' == ele['NameOfSubstance']:
                idf_ele_2 = idf_ele_2 + 1
                break
        for ele in ChemicalNameList:
            if 'Muscle, Skeletal' == ele['NameOfSubstance']:
                idf_ele_3 = idf_ele_3 + 1
                break
        for ele in ChemicalNameList:
            if 'Ryanodine Receptor Calcium Release Channel' == ele['NameOfSubstance']:
                idf_ele_4 = idf_ele_4 + 1
                break
        for ele in ChemicalNameList:
            if 'Calcium' == ele['NameOfSubstance']:
                idf_ele_5 = idf_ele_5 + 1
                break




        for eleM in MeshHeadingNameList:
            if 'Malignant Hyperthermia' == eleM['MeshHeadingName']:
                idf_eleM_1 = idf_eleM_1 + 1
                break
        for eleM in MeshHeadingNameList:
            if 'Muscle Contraction' == eleM['MeshHeadingName']:
                idf_eleM_2 = idf_eleM_2 + 1
                break
        for eleM in MeshHeadingNameList:
            if 'Muscle, Skeletal' == eleM['MeshHeadingName']:
                idf_eleM_3 = idf_eleM_3 + 1
                break
        for eleM in MeshHeadingNameList:
            if 'Ryanodine Receptor Calcium Release Channel' == eleM['MeshHeadingName']:
                idf_eleM_4 = idf_eleM_4 + 1
                break
        for eleM in MeshHeadingNameList:
            if re.findall(r'(Human|Humans)', eleM['MeshHeadingName']):
                idf_eleM_5 = idf_eleM_5 + 1
                break
        for eleM in MeshHeadingNameList:
            if 'Male' in eleM['MeshHeadingName']:
                idf_eleM_6 = idf_eleM_6 + 1
                break
        for eleM in MeshHeadingNameList:
            if 'Middle Aged' == eleM['MeshHeadingName']:
                idf_eleM_7 = idf_eleM_7 + 1
                break
        for eleM in MeshHeadingNameList:
            if 'Calcium' == eleM['MeshHeadingName']:
                idf_eleM_8 = idf_eleM_8 + 1
                break
        for eleM in MeshHeadingNameList:
            if re.findall(r'(Adult|Adults)', eleM['MeshHeadingName']):
                idf_eleM_9 = idf_eleM_9 + 1
                break

        for eleK in KeywordsList:
            if 'malignant hyperthermia' in str(eleK).lower():
                idf_eleK_1 = idf_eleK_1 + 1
                break
        for eleK in KeywordsList:
            if 'ryr1' in str(eleK).lower():
                idf_eleK_2 = idf_eleK_2 + 1
                break
        i=i+1
        print(i)

    print('idf_ele_1:'+str(idf_ele_1))
    print('idf_ele_2:'+str(idf_ele_2))
    print('idf_ele_3:'+str(idf_ele_3))
    print('idf_ele_4:'+str(idf_ele_4))
    print('idf_ele_5:' + str(idf_ele_5))
    print('idf_ele_6:'+str(idf_ele_6))
    print('idf_ele_7:' + str(idf_ele_7))
    print('idf_ele_8:' + str(idf_ele_8))

    print('idf_eleM_1:'+str(idf_eleM_1))
    print('idf_eleM_2:'+str(idf_eleM_2))
    print('idf_eleM_3:' + str(idf_eleM_3))
    print('idf_eleM_4:' + str(idf_eleM_4))
    print('idf_eleM_5:' + str(idf_eleM_5))
    print('idf_eleM_6:' + str(idf_eleM_6))
    print('idf_eleM_7:' + str(idf_eleM_7))
    print('idf_eleM_8:' + str(idf_eleM_8))
    print('idf_eleM_9:' + str(idf_eleM_9))
    print('idf_eleM_10:' + str(idf_eleM_10))
    print('idf_eleM_11:' + str(idf_eleM_11))
    print('idf_eleM_12:' + str(idf_eleM_12))

    print('idf_eleK_1:'+str(idf_eleK_1))
    print('idf_eleK_2:'+str(idf_eleK_2))
    print('idf_eleK_3:'+str(idf_eleK_3))


