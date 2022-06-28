
import pymongo

def caculate(ranked_list):
    sum=0
    total=0
    for j in ranked_list:
        if j in [1, 2]:
            sum = sum + 1
        if j in [1, 2,0]:
            total = total + 1
        if total==10:
            break
    return sum/total

myclient =pymongo.MongoClient("mongodb://localhost:27017/")
mytop=myclient['cs2017_result']


if __name__ == '__main__':
    sum_ng=0
    i=0
    for j in range(30):
        sorelist=[]
        temp1 = "cs2017_result_" + str(j + 1)
        mytop_100 = mytop[temp1]
        for x in mytop_100.find({}, {'related'}).limit(1000):
            sorelist.append(int(x['related']))
        sum = caculate(sorelist)
        sum_ng = sum_ng + sum
        i=i+1
        print(i)
    print(sum_ng/30)