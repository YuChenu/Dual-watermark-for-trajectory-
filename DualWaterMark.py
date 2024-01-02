import copy
import os
import pandas as pd
import tqdm
import math
from datetime import datetime
from shapely.geometry import geo
import numpy as np
from bitarray import bitarray
from tqdm import tqdm




def UniqueList(poly):
    uniquearr=[]
    for i in poly:
        if i not in uniquearr:
            uniquearr.append(i)
    return uniquearr

def point_to_line_Distance(point_a, point_b, point_c):
    
    if point_b[2] == point_c[2]:
        return -999
    else:
        slope = (point_b[3] - point_c[3]) / (point_b[2] - point_c[2])
        intercept = point_b[3] - slope * point_b[2]
        
        pointdistance = math.fabs(slope * point_a[2] - point_a[3] + intercept) / math.sqrt(1 + pow(slope, 2))
        return pointdistance

def sortFunction(date):
    return date[1]

def Time(time):
    return datetime.strptime(str(time), "%Y-%m-%d-%H:%M:%S")

def GetSEDDistance(x, y, x1, y1):
    return math.sqrt(math.pow((x - x1), 2) + math.pow((y - y1), 2))

def DouglasPeucker(Pointlist,threshold):

    resultlist=[]
    max_distance,max_distance_index=0,0
    
    for i in range(1,len(Pointlist)):
            point=Pointlist[i]
            distance = point_to_line_Distance(point, Pointlist[0], Pointlist[-1])
            if distance > max_distance:
                max_distance_index = i
                max_distance = distance


   
    if max_distance >= threshold:
        resultlist.append(Pointlist[max_distance_index])
        
        if max_distance_index-1!=0:
            resultleft=DouglasPeucker(Pointlist[0:max_distance_index],threshold)
            if len(resultleft)>0:
                resultlist.extend(resultleft)
        if len(Pointlist)-1!=max_distance_index:
            resultright=DouglasPeucker(Pointlist[max_distance_index:len(Pointlist)],threshold)
            if len(resultright)>0:
                resultlist.extend(resultright)
    else:
        if Pointlist[-1] not in resultlist:
            resultlist.append(Pointlist[-1])
        if Pointlist[0] not in resultlist:
            resultlist.append(Pointlist[0])
    return resultlist

def SimplificationbySEDTD2(tempGpoints, fromID, toID,threshold,result):
    
    maxindex = -1
    
    starttime = tempGpoints[fromID][1]
    endtime = tempGpoints[toID][1]
    t1 = datetime.strptime(str(starttime), "%Y-%m-%d-%H:%M:%S")
    t2 = datetime.strptime(str(endtime), "%Y-%m-%d-%H:%M:%S")
    totaltime = (t2 - t1).total_seconds()

    if toID-fromID>1 and totaltime>0:

        maxtimeratiodis = 0

        for i in range(fromID + 1, toID):
            ratio=(Time(tempGpoints[i][1]) - Time(tempGpoints[fromID][1])).total_seconds()/totaltime
            
            x_m = tempGpoints[fromID][2] + ratio * (tempGpoints[toID][2] - tempGpoints[fromID][2])
            y_m = tempGpoints[fromID][3] + ratio * (tempGpoints[toID][3] - tempGpoints[fromID][3])
            
            timeratiodis = GetSEDDistance(tempGpoints[i][2], tempGpoints[i][3], x_m, y_m)

            
            if timeratiodis >= maxtimeratiodis :
                maxtimeratiodis = timeratiodis
                maxindex = i

        
        if maxtimeratiodis >= threshold:
            result.append(tempGpoints[maxindex])
            SimplificationbySEDTD2(tempGpoints, fromID, maxindex, threshold,result)
            SimplificationbySEDTD2(tempGpoints, maxindex, toID, threshold,result)
        else:
            result.append(tempGpoints[fromID])
            result.append(tempGpoints[toID])
    else:
        return


def cal_ang(point_1, point_2, point_3):
    
    a = math.sqrt((point_2[0] - point_3[0]) * (point_2[0] - point_3[0]) + (point_2[1] - point_3[1]) * (
                point_2[1] - point_3[1]))
    b = math.sqrt((point_1[0] - point_3[0]) * (point_1[0] - point_3[0]) + (point_1[1] - point_3[1]) * (
                point_1[1] - point_3[1]))
    c = math.sqrt((point_1[0] - point_2[0]) * (point_1[0] - point_2[0]) + (point_1[1] - point_2[1]) * (
                point_1[1] - point_2[1]))
    if a == 0 or c == 0:
        print(point_1, point_2, point_3)
    B = math.degrees(math.acos((b * b - a * a - c * c) / (-2 * a * c)))
    return B



def GetCross(x1, y1, x2, y2, x, y):
    a = (x2 - x1, y2 - y1)
    b = (x - x1, y - y1)
    return a[0] * b[1] - a[1] * b[0]


def isInRectangle(x1, y1, x2, y2, x3, y3, x4, y4, x, y):
    return GetCross(x1, y1, x2, y2, x, y) * GetCross(x3, y3, x4, y4, x, y) >= 0 and GetCross(x2, y2, x3, y3, x,
                                                                                             y) * GetCross(x4,
                                                                                                           y4,
                                                                                                           x1,
                                                                                                           y1,
                                                                                                           x,
                                                                                                           y) >= 0



def GetFloat(data,lc):
    xs=str(data).split(".")[1]
    return int(xs[lc-1])



def str2bitarray(s):
    ret = bitarray(''.join([bin(int('1' + hex(c)[2:], 16))[3:] for c in s.encode('utf-8')]))
    return ret



def bitarray2str(bit):
    return bit.tobytes().decode('utf-8')


def cal_disratio(point1, point2, point3):
   
    dis1 = math.sqrt(math.pow(point1[0] - point2[0], 2) + math.pow(point1[1] - point2[1], 2))
    dis2 = math.sqrt(math.pow(point3[0] - point2[0], 2) + math.pow(point3[1] - point2[1], 2))
    return dis1 / dis2



def QIM(data, N, bit, lc):
    
    xiaoshu = str(data).split(".")[1]
    single_bit = xiaoshu[lc - 1]

    new_single_bit = int(int(single_bit) / N) * N + bit

    xiaoshu = xiaoshu[:lc - 1] + str(new_single_bit) + xiaoshu[lc:len(xiaoshu)]

    new_data = int(data) + int(xiaoshu) / math.pow(10, len(xiaoshu))
    return new_data



def RQIM(data,N):
    return data-(int(int(data) / N))*N



def GetFloatlocBit(data, lc):
    xiaoshu = str(data).split(".")[1][lc - 1]
    return xiaoshu



def NC(arr1, arr2):
    Sum = 0
    suma1 = 0
    suma2 = 0
    for i in range(len(arr1)):
        Sum = Sum + (arr1[i] * arr2[i])
        suma1 = suma1 + (arr1[i] * arr1[i])
        suma2 = suma2 + (arr2[i] * arr2[i])
    cor = Sum / (math.sqrt(suma1) * math.sqrt(suma2))
    return cor



def RoatePoint(p1, A, B, oldangle, newangle):
    
    i = B[0] - A[0]
    j = B[1] - A[1]

    d_p1p2 = math.sqrt((A[0] - p1[0]) * (A[0] - p1[0]) + (A[1] - p1[1]) * (A[1] - p1[1]))
    d_p2p3 = math.sqrt((A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1] - B[1]))

    t = (math.cos(math.radians(newangle)) * d_p1p2 * d_p2p3 + (i * A[0])) / j
    a = 1 + (i * i) / (j * j)
    b = -(2 * A[0] + 2 * t * (i / j))
    c = math.pow(A[0], 2) + math.pow(t, 2) - math.pow(d_p1p2, 2)

    x1 = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
    y1 = t + A[1] - (i / j * x1)

    x2 = (-b - math.sqrt(b * b - 4 * a * c)) / (2 * a)
    y2 = t + A[1] - (i / j * x2)

    cor_x_len = len(str(p1[0]).split(".")[1])
    cor_y_len = len(str(p1[1]).split(".")[1])

    # print(x1,y1,round(cal_ang((x1,y1),A,B),len(str(newangle).split(".")[1])))
    # print(x2,y2,round(cal_ang((x2,y2),A,B),len(str(newangle).split(".")[1])))
    # print(x1,y1,round(cal_ang((round(x1,cor_x_len),round(y1,cor_y_len)),A,B),len(str(newangle).split(".")[1])))
    # print(x2,y2,round(cal_ang((round(x2,cor_x_len),round(y2,cor_y_len)),A,B),len(str(newangle).split(".")[1])))

    
    if newangle == oldangle:
        return p1
    else:
        CAC1_1 = cal_ang((round(x1, cor_x_len), round(y1, cor_y_len)), A, p1)
        CAC1_2 = cal_ang((round(x2, cor_x_len), round(y2, cor_y_len)), A, p1)
        if CAC1_1 < CAC1_2:
            return (x1, y1)
        if CAC1_1 >= CAC1_2:
            return (x2, y2)



def datetimeToint(datetime):
    year = int(str(datetime).split("-")[0])
    month = int(str(datetime).split("-")[1])
    day = int(str(datetime).split("-")[2])
    time = str(datetime).split("-")[3]
    hour = int(time.split(":")[0])
    minute = int(time.split(":")[1])
    second = int(time.split(":")[2])

    datetimelist = [year, month, day, hour, minute, second]
    sum = 0
    for t in datetimelist:
        sum += 5 + 11 * t
    return int(sum / 2) % 2



def ZipPoint(p1, A, B, k):
    slope_AC = (A[1] - p1[1]) / (A[0] - p1[0])
    dis_AB = math.sqrt(math.pow(A[0] - B[0], 2) + math.pow(A[1] - B[1], 2))
    a = 1 + slope_AC * slope_AC
    b = -(2 * slope_AC * slope_AC * A[0] + 2 * A[0])
    c = (1 + slope_AC * slope_AC) * (A[0] * A[0]) - (
            k * k * (math.pow(A[0] - B[0], 2) + math.pow(A[1] - B[1], 2)))

    delt = b * b - (4 * a * c)
    x1 = (-b + math.sqrt(delt)) / (2 * a)
    y1 = slope_AC * x1 + A[1] - slope_AC * A[0]

    x2 = (-b - math.sqrt(delt)) / (2 * a)
    y2 = slope_AC * x2 + A[1] - slope_AC * A[0]

    dis_C1p1 = math.sqrt(math.pow(x1 - p1[0], 2) + math.pow(y1 - p1[1], 2))
    dis_C1p2 = math.sqrt(math.pow(x2 - p1[0], 2) + math.pow(y2 - p1[1], 2))

    if dis_C1p1 < dis_C1p2:
        # return (round(x1,6),round(y1,6))
        return (x1, y1)
    if dis_C1p2 < dis_C1p1:
        # return (round(x2,6),round(y2,6))
        return (x2, y2)



def Selection_sort_bydis(A, firstpoint, secondpoint):
    for i in range(len(A)):
        min_idx = i
        for j in range(i + 1, len(A)):
            if point_to_line_Distance(A[min_idx], firstpoint, secondpoint) > point_to_line_Distance(A[j],
                                                                                                    firstpoint,
                                                                                                    secondpoint):
                min_idx = j
        A[i], A[min_idx] = A[min_idx], A[i]


# Embed Roubust WaterMark
def RoubustWaterMarkEmbed(data, bit, QIMStep, QIMlc,firstmindis_point,secondmindis_point):
    locationlist = []
    for count in range(len(data)):
        point = data[count]
        
        dis_ratio = cal_disratio((point[2], point[3]), (firstmindis_point[2], firstmindis_point[3]),
                                 (secondmindis_point[2], secondmindis_point[3]))

        
        if point[4] == "T":
            
            angel_fp = cal_ang((point[2], point[3]), (firstmindis_point[2], firstmindis_point[3]),
                               (secondmindis_point[2], secondmindis_point[3]))

            
            listindex = list(str(int(angel_fp)) + str(int(dis_ratio)))
            a = 9
            b = 11
            sum = 0
            for k in listindex:
                sum = sum + (int(k) * a + b)

            location = ((int(angel_fp) + int(dis_ratio)) * a + sum * b) % len(bit)

            locationlist.append(location)

            new_angle_QIM = QIM(angel_fp, QIMStep, bit[location], QIMlc)

           
            new_point_cor = RoatePoint((point[2], point[3]), (firstmindis_point[2], firstmindis_point[3]),
                                       (secondmindis_point[2], secondmindis_point[3]), angel_fp, new_angle_QIM)

            new_angle_embed = cal_ang(new_point_cor, (firstmindis_point[2], firstmindis_point[3]),
                                      (secondmindis_point[2], secondmindis_point[3]))

            if (GetFloatlocBit(new_angle_QIM, QIMlc) == GetFloatlocBit(new_angle_embed, QIMlc)) == False:
                print(point, "鲁棒水印嵌入失败")
            point[2] = new_point_cor[0]
            point[3] = new_point_cor[1]
    
    embed_watermark = [0] * len(bit)
    for i in locationlist:
        embed_watermark[i] = bit[i]
    
    return NC(bit, embed_watermark)



def FragileWaterMarkEmbed(data, newTraID, QIMStep, QIMlc,firstmindis_point,secondmindis_point):

    for count in range(len(data)):
        
        if count != 0:
            point = data[count]
            beforepoint = data[count - 1]

            
            dateinfor = datetimeToint(beforepoint[1])

            
            dis_ratio = cal_disratio((point[2], point[3]), (firstmindis_point[2], firstmindis_point[3]),
                                     (secondmindis_point[2], secondmindis_point[3]))

            
            new_dis_ratio_QIM = QIM(dis_ratio, QIMStep, dateinfor, QIMlc)
            
            new_point = ZipPoint((point[2], point[3]), (firstmindis_point[2], firstmindis_point[3]),
                                 (secondmindis_point[2], secondmindis_point[3]), new_dis_ratio_QIM)

            
            new_dis_ratio_embed = cal_disratio(new_point, (firstmindis_point[2], firstmindis_point[3]),
                                               (secondmindis_point[2], secondmindis_point[3]))

            if (GetFloatlocBit(new_dis_ratio_QIM, QIMlc) == GetFloatlocBit(new_dis_ratio_embed,
                                                                           QIMlc)) == False:
                print(point, new_point, dis_ratio, new_dis_ratio_QIM, new_dis_ratio_embed, "脆弱水印嵌入失败")
            point[2] = new_point[0]
            point[3] = new_point[1]





def sortFunction(date):
    return date[1]




def RoubustWaterMarkExtract(data, emded_bit, QIMStep, QIMlc,firstmindis_point,secondmindis_point):

    
    def MinMajorVote(num):
        major = num[0]
        count = 0
        for i in num:
            if count == 0:
                major = i
            if i == major:
                count += 1
            else:
                count -= 1
        if num.count(major) > len(num) / 2:
            return major
        else:
           
            return 1

    
    def NC(arr1, arr2):
        Sum = 0
        suma1 = 0
        suma2 = 0
        for i in range(len(arr1)):
            if arr1[i] != -1:
                Sum = Sum + (arr1[i] * arr2[i])
                suma1 = suma1 + (arr1[i] * arr1[i])
                suma2 = suma2 + (arr2[i] * arr2[i])
        cor = Sum / (math.sqrt(suma1) * math.sqrt(suma2))
        return cor

    waterbitlocation = [[] for i in range(len(emded_bit))]

    for count in range(len(data)):
        point = data[count]
        
        dis_ratio = cal_disratio((point[2], point[3]), (firstmindis_point[2], firstmindis_point[3]),
                                 (secondmindis_point[2], secondmindis_point[3]))

        
        if point[4] == "T":
            
            angel_fp = cal_ang((point[2], point[3]), (firstmindis_point[2], firstmindis_point[3]),
                               (secondmindis_point[2], secondmindis_point[3]))

           
            listindex = list(str(int(angel_fp)) + str(int(dis_ratio)))
            a = 9
            b = 11
            sum = 0
            for k in listindex:
                sum = sum + (int(k) * a + b)

            location = ((int(angel_fp) + int(dis_ratio)) * a + sum * b) % len(emded_bit)



            
            new_point_cor = cal_ang((point[2], point[3]), (firstmindis_point[2], firstmindis_point[3]),
                                       (secondmindis_point[2], secondmindis_point[3]))

            extract_single_bit = RQIM(GetFloat(new_point_cor,QIMlc), QIMStep)
            waterbitlocation[location].append(extract_single_bit)


    extract_bit=[]
    for i in range(len(waterbitlocation)):
        if len(waterbitlocation[i]) > 0:
            majorbit = MinMajorVote(waterbitlocation[i])
            extract_bit.append(majorbit)
        else:
            majorbit = -1
            extract_bit.append(majorbit)




    # print(locationlist,embed_bit,embed_watermark)
    return NC(emded_bit, extract_bit)




def FragileWaterMarkExtract(data, newTraID, QIMStep, QIMlc,firstmindis_point,secondmindis_point):
   

    extract_result=[]
    for count in range(len(data)):
        
        if count != 0:
            point = data[count]
            beforepoint = data[count - 1]

            
            dateinfor = datetimeToint(beforepoint[1])

            # 计算距离比值
            dis_ratio = cal_disratio((point[2], point[3]), (firstmindis_point[2], firstmindis_point[3]),
                                     (secondmindis_point[2], secondmindis_point[3]))

            
            new_dis_ratio_RQIM = RQIM(GetFloat(dis_ratio,QIMlc), QIMStep)

            
            if dateinfor==new_dis_ratio_RQIM:
                extract_result.append([data[count],True])
            else:
                extract_result.append([data[count],False])

    return extract_result







def WaterMarkEmbed(tra_all,dpthreshold,tdtrthreshold,robust_threshold,fragile_threshold,embedbit):

    
    tra_new = []
    tra_tmp = []
    for single_tra in tqdm(tra_all):
        for i in range(len(single_tra)):
            if i < len(single_tra) - 1 and abs(
                    datetime.strptime(single_tra[i + 1][1], "%Y-%m-%d-%H:%M:%S").second - datetime.strptime(single_tra[i][1],
                                                                                                         "%Y-%m-%d-%H:%M:%S").second) < 60:
                tra_tmp.append(single_tra[i])
            else:
                tra_tmp.append(single_tra[i])
                tra_new.append(tra_tmp.copy())
                tra_tmp.clear()



   


    Featurepoint = []
    for sub_tra in tra_new:

       
        DPPointZip = DouglasPeucker(sub_tra, dpthreshold)

       
        sub_tra.sort(key=sortFunction)
        TDTRPointZip = []
        SimplificationbySEDTD2(sub_tra, 0, len(sub_tra) - 1, tdtrthreshold, TDTRPointZip)

        
        UniqueDPPoint = UniqueList(DPPointZip)
        UniqueTDTRPoint = UniqueList(TDTRPointZip)

        if len(UniqueDPPoint) <= len(UniqueTDTRPoint):
            Featurepointtmp = [v for v in UniqueDPPoint if v in UniqueTDTRPoint]
        else:
            Featurepointtmp = [v for v in UniqueTDTRPoint if v in UniqueDPPoint]

        for i in Featurepointtmp:
            Featurepoint.append(i)

    
    for i in tra_all:
        if i in Featurepoint:
            i[4] = "T"


    
    Featurepoint_onlyCor = []
    for i in Featurepoint:
        Featurepoint_onlyCor.append([i[2], i[3]])


   
    if len(Featurepoint_onlyCor) > 2:
        multipoint = geo.MultiPoint(np.array(Featurepoint_onlyCor))
        min_rect = multipoint.minimum_rotated_rectangle
        pointlist = str(min_rect)[10:-2].split(",")

        min_rect_r = []
        for i in range(len(pointlist) - 1):

            if pointlist[i][:1] == " ":
                cor = pointlist[i][1:].split(" ")
            else:
                cor = pointlist[i].split(" ")
            min_rect_r.append([float(cor[0]), float(cor[1])])

   
        maxdiagonaldis = 0  
        maxindex = -1  
        maxpoint = []  
        for point_num in range(len(min_rect_r)):
            dis = math.sqrt(math.pow(min_rect_r[0][0] - min_rect_r[point_num][0], 2) + math.pow(
                min_rect_r[0][1] - min_rect_r[point_num][1], 2))
            if dis > maxdiagonaldis:
                maxdiagonaldis = dis
                maxpoint = min_rect_r[point_num]
                maxindex = point_num

       
        centerpoint = [(min_rect_r[0][0] + maxpoint[0]) / 2, (min_rect_r[0][1] + maxpoint[1]) / 2]

       
        firstmindis, secondmindis, firstmindis_point, secondmindis_point = 99999, 199999, [0, 0], [0, 0]
        for feap in Featurepoint:
            feap_x = feap[2]
            feap_y = feap[3]
            dis = math.sqrt(math.pow(float(feap_x) - centerpoint[0], 2) + math.pow(float(feap_y) - centerpoint[1], 2))

            if dis < firstmindis:
                secondmindis = firstmindis
                secondmindis_point = firstmindis_point
                firstmindis = dis
                firstmindis_point = feap
            else:
                if dis < secondmindis:
                    secondmindis = dis
                    secondmindis_point = feap


       


       
        right_tra, left_tra = [], []
        for i in tra_all:
            if i[2] != firstmindis_point[2] and i[2] != secondmindis_point[2] and i[3] != firstmindis_point[3] and i[3] != \
                    secondmindis_point[3]:
                ratio = cal_ang((i[2], i[3]), (firstmindis_point[2], firstmindis_point[3]),
                                (secondmindis_point[2], secondmindis_point[3]))
                if 0 < ratio < 90:
                    right_tra.append(i)
                if 90 < ratio < 180:
                    left_tra.append(i)

       

        
        def GetRightTra_Vertex(tra, pointA_x, pointA_y, pointB_x, pointB_y):
            
            maxdis_hori = 0
            maxdis_ver = 0
            slope_AB = (pointB_y - pointA_y) / (pointB_x - pointA_x)
            slope_AC = -(1 / slope_AB)
            b_AC = pointA_y - slope_AC * pointA_x

            for right_tra_p in tra:
                
                p_to_AB = point_to_line_Distance(right_tra_p, (0, 0, pointA_x, pointA_y, 0),
                                                 (0, 0, pointB_x, pointB_y, 0))
               
                p_to_AC = math.fabs(
                    (slope_AC * right_tra_p[2] - right_tra_p[3] + b_AC) / (math.sqrt(1 + math.pow(slope_AC, 2))))

                
                if p_to_AB > maxdis_hori:
                    maxdis_hori = p_to_AB
                if p_to_AC > maxdis_ver:
                    maxdis_ver = p_to_AC

            if maxdis_ver > maxdis_hori:
                square_len = maxdis_ver + 0.001
            else:
                square_len = maxdis_hori + 0.001

            
            delt_a = 1 + math.pow(slope_AC, 2)

            delt_b = -(2 * pointA_x - 2 * slope_AC * (b_AC - pointA_y))

            delt_c = math.pow(b_AC - pointA_y, 2) - math.pow(square_len, 2) + math.pow(pointA_x, 2)

            delt = math.sqrt(delt_b * delt_b - (4 * delt_a * delt_c))

        
            C1_x = (-delt_b + delt) / (2 * delt_a)
            C1_y = slope_AC * C1_x + b_AC

            C2_x = (-delt_b - math.sqrt(delt_b * delt_b - 4 * delt_a * delt_c)) / (2 * delt_a)
            C2_y = slope_AC * C2_x + b_AC

            

            def Get_parallelC_pointRight(slope_AC, slope_AB, C_X, C_Y, dissquare, firstmindis_point_x,
                                         firstmindis_point_y, secondmindis_point_x, secondmindis_point_y):
                
                t1_C1 = slope_AC - slope_AB
                t2_C1 = (slope_AB - slope_AC) * C_X
                t3_C1 = dissquare * dissquare * (1 + slope_AC * slope_AC)
                x1_C2 = ((-(2 * t1_C1 * t2_C1)) + math.sqrt(
                    4 * t1_C1 * t1_C1 * t2_C1 * t2_C1 - 4 * t1_C1 * t1_C1 * (t2_C1 * t2_C1 - t3_C1))) / (
                                    2 * t1_C1 * t1_C1)
                y1_C2 = slope_AB * x1_C2 + C_Y - slope_AB * C_X

                x2_C2 = ((-(2 * t1_C1 * t2_C1)) - math.sqrt(
                    4 * t1_C1 * t1_C1 * t2_C1 * t2_C1 - 4 * t1_C1 * t1_C1 * (t2_C1 * t2_C1 - t3_C1))) / (
                                    2 * t1_C1 * t1_C1)
                y2_C2 = slope_AB * x2_C2 + C_Y - slope_AB * C_X

                ratio_x1_C2 = cal_ang((x1_C2, y1_C2), (firstmindis_point_x, firstmindis_point_y),
                                      (secondmindis_point_x, secondmindis_point_y))
                ratio_x2_C2 = cal_ang((x2_C2, y2_C2), (firstmindis_point_x, firstmindis_point_y),
                                      (secondmindis_point_x, secondmindis_point_y))
                if ratio_x1_C2 < 90:
                    return x1_C2, y1_C2
                if ratio_x2_C2 < 90:
                    return x2_C2, y2_C2

            D1_x, D1_y = Get_parallelC_pointRight(slope_AC, slope_AB, C1_x, C1_y, square_len, pointA_x, pointA_y,
                                                  pointB_x, pointB_y)
            D2_x, D2_y = Get_parallelC_pointRight(slope_AC, slope_AB, C2_x, C2_y, square_len, pointA_x, pointA_y,
                                                  pointB_x, pointB_y)

            return [C1_x, C1_y], [D1_x, D1_y], [D2_x, D2_y], [C2_x, C2_y]

        def GetLeftTra_Vertex(tra, pointA_x, pointA_y, pointB_x, pointB_y):
            
            maxdis_hori = 0
            maxdis_ver = 0
            slope_AB = (pointB_y - pointA_y) / (pointB_x - pointA_x)
            slope_AC = -(1 / slope_AB)
            b_AC = pointA_y - slope_AC * pointA_x

            for right_tra_p in tra:
                
                p_to_AB = point_to_line_Distance(right_tra_p, (0, 0, pointA_x, pointA_y, 0),
                                                 (0, 0, pointB_x, pointB_y, 0))
                
                p_to_AC = math.fabs(
                    (slope_AC * right_tra_p[2] - right_tra_p[3] + b_AC) / (math.sqrt(1 + math.pow(slope_AC, 2))))

                
                if p_to_AB > maxdis_hori:
                    maxdis_hori = p_to_AB
                if p_to_AC > maxdis_ver:
                    maxdis_ver = p_to_AC

            if maxdis_ver > maxdis_hori:
                square_len = maxdis_ver
            else:
                square_len = maxdis_hori

            
            delt_a = 1 + math.pow(slope_AC, 2)

            delt_b = -(2 * pointA_x - 2 * slope_AC * (b_AC - pointA_y))

            delt_c = math.pow(b_AC - pointA_y, 2) - math.pow(square_len, 2) + math.pow(pointA_x, 2)

            delt = math.sqrt(delt_b * delt_b - (4 * delt_a * delt_c))

           
            C1_x = (-delt_b + delt) / (2 * delt_a)
            C1_y = slope_AC * C1_x + b_AC

            C2_x = (-delt_b - math.sqrt(delt_b * delt_b - 4 * delt_a * delt_c)) / (2 * delt_a)
            C2_y = slope_AC * C2_x + b_AC

            

            def Get_parallelC_pointLeft(slope_AC, slope_AB, C_X, C_Y, dissquare, firstmindis_point_x,
                                        firstmindis_point_y, secondmindis_point_x, secondmindis_point_y):
                
                t1_C1 = slope_AC - slope_AB
                t2_C1 = (slope_AB - slope_AC) * C_X
                t3_C1 = dissquare * dissquare * (1 + slope_AC * slope_AC)
                x1_C2 = ((-(2 * t1_C1 * t2_C1)) + math.sqrt(
                    4 * t1_C1 * t1_C1 * t2_C1 * t2_C1 - 4 * t1_C1 * t1_C1 * (t2_C1 * t2_C1 - t3_C1))) / (
                                    2 * t1_C1 * t1_C1)
                y1_C2 = slope_AB * x1_C2 + C_Y - slope_AB * C_X

                x2_C2 = ((-(2 * t1_C1 * t2_C1)) - math.sqrt(
                    4 * t1_C1 * t1_C1 * t2_C1 * t2_C1 - 4 * t1_C1 * t1_C1 * (t2_C1 * t2_C1 - t3_C1))) / (
                                    2 * t1_C1 * t1_C1)
                y2_C2 = slope_AB * x2_C2 + C_Y - slope_AB * C_X

                ratio_x1_C2 = cal_ang((x1_C2, y1_C2), (firstmindis_point_x, firstmindis_point_y),
                                      (secondmindis_point_x, secondmindis_point_y))
                ratio_x2_C2 = cal_ang((x2_C2, y2_C2), (firstmindis_point_x, firstmindis_point_y),
                                      (secondmindis_point_x, secondmindis_point_y))
                if ratio_x1_C2 > 90:
                    return x1_C2, y1_C2
                if ratio_x2_C2 > 90:
                    return x2_C2, y2_C2

            D1_x, D1_y = Get_parallelC_pointLeft(slope_AC, slope_AB, C1_x, C1_y, square_len, pointA_x, pointA_y,
                                                 pointB_x, pointB_y)
            D2_x, D2_y = Get_parallelC_pointLeft(slope_AC, slope_AB, C2_x, C2_y, square_len, pointA_x, pointA_y,
                                                 pointB_x, pointB_y)
            return [C1_x, C1_y], [D1_x, D1_y], [D2_x, D2_y], [C2_x, C2_y]

        rightvec1, rightvec2, rightvec3, rightvec4 = GetRightTra_Vertex(right_tra, firstmindis_point[2],
                                                                        firstmindis_point[3], secondmindis_point[2],
                                                                        secondmindis_point[3])
        leftvec1, leftvec2, leftvec3, leftvec4 = GetLeftTra_Vertex(left_tra, firstmindis_point[2], firstmindis_point[3],
                                                                   secondmindis_point[2], secondmindis_point[3])

        vec_right = [rightvec1, rightvec2, rightvec3, rightvec4]
        vec_left = [leftvec1, leftvec2, leftvec3, leftvec4]

        
        zone_fp=[]

        
        def DivideRecRobust(veccor, point, flag, threshold, robuest_zone):
            
            inAreaPoint = []
            for p in point:
                if isInRectangle(veccor[0][0], veccor[0][1], veccor[1][0], veccor[1][1], veccor[2][0], veccor[2][1],
                                 veccor[3][0], veccor[3][1], p[2], p[3]):
                    inAreaPoint.append(p)

            
            featurepoint_count = 0
            featurepointlist = []
            for feature_point in inAreaPoint:
                
                if feature_point[4] == "T":
                    featurepoint_count += 1
                    featurepointlist.append(feature_point)

           
            if featurepoint_count >= threshold:
                
                if flag == True:
                    newveccor1 = [veccor[0], veccor[1],
                                  [(veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2],
                                  [(veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2]]
                    newveccor2 = [[(veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2],
                                  [(veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2], veccor[2],
                                  veccor[3]]
                    newveccor1count, newveccor2count = 0, 0
                    for p in featurepointlist:
                        if isInRectangle(veccor[0][0], veccor[0][1], veccor[1][0], veccor[1][1],
                                         (veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2,
                                         (veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2, p[2],
                                         p[3]):
                            newveccor1count += 1
                        if isInRectangle((veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2,
                                         (veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2,
                                         veccor[2][0], veccor[2][1], veccor[3][0], veccor[3][1], p[2], p[3]):
                            newveccor2count += 1
                   
                    if newveccor1count < threshold and newveccor2count < threshold:
                        robuest_zone.append([inAreaPoint, veccor])
                    else:
                        DivideRecRobust(newveccor1, inAreaPoint, False, threshold, robuest_zone)
                        DivideRecRobust(newveccor2, inAreaPoint, False, threshold, robuest_zone)

                if flag == False:
                    newveccor1 = [veccor[0], [(veccor[0][0] + veccor[1][0]) / 2, (veccor[0][1] + veccor[1][1]) / 2],
                                  [(veccor[2][0] + veccor[3][0]) / 2, (veccor[2][1] + veccor[3][1]) / 2], veccor[3]]
                    newveccor2 = [[(veccor[0][0] + veccor[1][0]) / 2, (veccor[0][1] + veccor[1][1]) / 2], veccor[1],
                                  veccor[2], [(veccor[2][0] + veccor[3][0]) / 2, (veccor[2][1] + veccor[3][1]) / 2]]
                    newveccor1count, newveccor2count = 0, 0
                    for p in featurepointlist:
                        if isInRectangle(veccor[0][0], veccor[0][1], (veccor[0][0] + veccor[1][0]) / 2,
                                         (veccor[0][1] + veccor[1][1]) / 2, (veccor[2][0] + veccor[3][0]) / 2,
                                         (veccor[2][1] + veccor[3][1]) / 2, veccor[3][0], veccor[3][1], p[2], p[3]):
                            newveccor1count += 1
                        if isInRectangle((veccor[0][0] + veccor[1][0]) / 2, (veccor[0][1] + veccor[1][1]) / 2,
                                         veccor[1][0], veccor[1][1], veccor[2][0], veccor[2][1],
                                         (veccor[2][0] + veccor[3][0]) / 2, (veccor[2][1] + veccor[3][1]) / 2, p[2],
                                         p[3]):
                            newveccor2count += 1
                   
                    if newveccor1count < threshold and newveccor2count < threshold:
                        robuest_zone.append([inAreaPoint, veccor])
                    else:
                        DivideRecRobust(newveccor1, inAreaPoint, True, threshold, robuest_zone)
                        DivideRecRobust(newveccor2, inAreaPoint, True, threshold, robuest_zone)
            else:
                if len(inAreaPoint) > 0:
                    robuest_zone.append([inAreaPoint, veccor])

        def area_of_quadrangle(x1, y1, x2, y2, x3, y3, x4, y4):
            
            area1 = abs(0.5 * (x1 * y2 + x2 * y3 + x3 * y4 + x4 * y1 - x2 * y1 - x3 * y2 - x4 * y3 - x1 * y4))
            area2 = abs(0.5 * (x1 * y3 + x3 * y4 + x4 * y2 + x2 * y1 - x3 * y1 - x4 * y2 - x2 * y3 - x1 * y4))
            
            area = area1 + area2
            return area

        def DivedeRecFragile(veccor, point, flag, threshold, zone):
            inAreaPoint = []
            for p in point:
                if isInRectangle(veccor[0][0], veccor[0][1], veccor[1][0], veccor[1][1], veccor[2][0], veccor[2][1],
                                 veccor[3][0], veccor[3][1], p[2], p[3]):
                    inAreaPoint.append(p)

            
            recarea = area_of_quadrangle(veccor[0][0], veccor[0][1], veccor[1][0], veccor[1][1], veccor[2][0],
                                         veccor[2][1], veccor[3][0], veccor[3][1])
            # print(recarea)

            
            if recarea >= threshold and len(inAreaPoint) > 0:
               
                if flag == True:
                    newveccor1 = [veccor[0], veccor[1],
                                  [(veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2],
                                  [(veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2]]
                    newveccor2 = [[(veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2],
                                  [(veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2], veccor[2],
                                  veccor[3]]

                    a1 = area_of_quadrangle(veccor[0][0], veccor[0][1], veccor[1][0], veccor[1][1],
                                            (veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2,
                                            (veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2)

                    a2 = area_of_quadrangle((veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2,
                                            (veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2,
                                            veccor[2][0], veccor[2][1], veccor[3][0], veccor[3][1])

                    
                    if a1 < threshold and a2 < threshold:
                        zone.append(inAreaPoint)
                    else:
                        DivedeRecFragile(newveccor1, inAreaPoint, False, threshold, zone)
                        DivedeRecFragile(newveccor2, inAreaPoint, False, threshold, zone)

                if flag == False:
                    newveccor1 = [veccor[0], [(veccor[0][0] + veccor[1][0]) / 2, (veccor[0][1] + veccor[1][1]) / 2],
                                  [(veccor[2][0] + veccor[3][0]) / 2, (veccor[2][1] + veccor[3][1]) / 2], veccor[3]]
                    newveccor2 = [[(veccor[0][0] + veccor[1][0]) / 2, (veccor[0][1] + veccor[1][1]) / 2], veccor[1],
                                  veccor[2], [(veccor[2][0] + veccor[3][0]) / 2, (veccor[2][1] + veccor[3][1]) / 2]]

                    a3 = area_of_quadrangle(veccor[0][0], veccor[0][1], (veccor[0][0] + veccor[1][0]) / 2,
                                            (veccor[0][1] + veccor[1][1]) / 2, (veccor[2][0] + veccor[3][0]) / 2,
                                            (veccor[2][1] + veccor[3][1]) / 2, veccor[3][0], veccor[3][1])

                    a4 = area_of_quadrangle((veccor[0][0] + veccor[1][0]) / 2, (veccor[0][1] + veccor[1][1]) / 2,
                                            veccor[1][0], veccor[1][1], veccor[2][0], veccor[2][1],
                                            (veccor[2][0] + veccor[3][0]) / 2, (veccor[2][1] + veccor[3][1]) / 2)

                    
                    if a3 < threshold and a4 < threshold:
                        zone.append(inAreaPoint)
                    else:
                        DivedeRecFragile(newveccor1, inAreaPoint, True, threshold, zone)
                        DivedeRecFragile(newveccor2, inAreaPoint, True, threshold, zone)
            else:
                if len(inAreaPoint) > 0:
                    zone.append(inAreaPoint)

       
        DivideRecRobust(vec_right, right_tra, True, robust_threshold,zone_fp)
        DivideRecRobust(vec_left, left_tra, True, robust_threshold,zone_fp)

        final_zone = []
        areathreshold = 0.00001

        for z in zone_fp:
            zone_tmp = z[0]

            count = 0
            for i in zone_tmp:
                if i[4] == "T":
                    count = count + 1
            # print(count)
            if robust_threshold < count:
                DivedeRecFragile(z[1], zone_tmp, True, fragile_threshold, final_zone)
            else:
                final_zone.append(zone_tmp)


        

        
        embed_bit = str2bitarray(embedbit).tolist()

        slope_AB = (secondmindis_point[3] - firstmindis_point[3]) / (secondmindis_point[2] - firstmindis_point[2])
        b_AB = firstmindis_point[3] - (slope_AB * firstmindis_point[2])

       
        new_tra = []

        for zone in final_zone:
            changedata = copy.deepcopy(zone)

            
            Selection_sort_bydis(changedata, firstmindis_point, secondmindis_point)


            fpcount = 0
            for fp in zone:
                if fp[4] == "T":
                    fpcount += 1

           
            if fpcount >= robust_threshold:
                RoubustWaterMarkEmbed(changedata, embed_bit, 2, 1,firstmindis_point,secondmindis_point)

          
            FragileWaterMarkEmbed(changedata, "ID", 2, 1,firstmindis_point,secondmindis_point)

            for t in changedata:
                new_tra.append(t)

       
        new_tra.append(firstmindis_point)
        new_tra.append(secondmindis_point)

        
        new_tra.sort(key=sortFunction)

        return new_tra



def WaterMarkExtract(tra_all,dpthreshold,tdtrthreshold,robust_threshold,fragile_threshold,embedbit):
    
    tra_new = []
    tra_tmp = []
    for i in range(len(tra_all)):
        if i < len(tra_all) - 1 and abs(
                datetime.strptime(tra_all[i + 1][1], "%Y-%m-%d-%H:%M:%S").second - datetime.strptime(tra_all[i][1],
                                                                                                     "%Y-%m-%d-%H:%M:%S").second) < 60:
            tra_tmp.append(tra_all[i])
        else:
            tra_tmp.append(tra_all[i])
            tra_new.append(tra_tmp.copy())
            tra_tmp.clear()



    Featurepoint = []
    for sub_tra in tra_new:

        
        DPPointZip = DouglasPeucker(sub_tra, dpthreshold)

        
        sub_tra.sort(key=sortFunction)
        TDTRPointZip = []
        SimplificationbySEDTD2(sub_tra, 0, len(sub_tra) - 1, tdtrthreshold, TDTRPointZip)

       
        UniqueDPPoint = UniqueList(DPPointZip)
        UniqueTDTRPoint = UniqueList(TDTRPointZip)

        if len(UniqueDPPoint) <= len(UniqueTDTRPoint):
            Featurepointtmp = [v for v in UniqueDPPoint if v in UniqueTDTRPoint]
        else:
            Featurepointtmp = [v for v in UniqueTDTRPoint if v in UniqueDPPoint]

        for i in Featurepointtmp:
            Featurepoint.append(i)

    
    for i in tra_all:
        if i in Featurepoint:
            i[4] = "T"


    
    Featurepoint_onlyCor = []
    for i in Featurepoint:
        Featurepoint_onlyCor.append([i[2], i[3]])


    
    if len(Featurepoint_onlyCor) > 2:
        multipoint = geo.MultiPoint(np.array(Featurepoint_onlyCor))
        min_rect = multipoint.minimum_rotated_rectangle
        pointlist = str(min_rect)[10:-2].split(",")

        min_rect_r = []
        for i in range(len(pointlist) - 1):

            if pointlist[i][:1] == " ":
                cor = pointlist[i][1:].split(" ")
            else:
                cor = pointlist[i].split(" ")
            min_rect_r.append([float(cor[0]), float(cor[1])])



        
        maxdiagonaldis = 0  
        maxindex = -1  
        maxpoint = []  
        for point_num in range(len(min_rect_r)):
            dis = math.sqrt(math.pow(min_rect_r[0][0] - min_rect_r[point_num][0], 2) + math.pow(
                min_rect_r[0][1] - min_rect_r[point_num][1], 2))
            if dis > maxdiagonaldis:
                maxdiagonaldis = dis
                maxpoint = min_rect_r[point_num]
                maxindex = point_num

        
        centerpoint = [(min_rect_r[0][0] + maxpoint[0]) / 2, (min_rect_r[0][1] + maxpoint[1]) / 2]

        
        firstmindis, secondmindis, firstmindis_point, secondmindis_point = 99999, 199999, [0, 0], [0, 0]
        for feap in Featurepoint:
            feap_x = feap[2]
            feap_y = feap[3]
            dis = math.sqrt(math.pow(float(feap_x) - centerpoint[0], 2) + math.pow(float(feap_y) - centerpoint[1], 2))

            if dis < firstmindis:
                secondmindis = firstmindis
                secondmindis_point = firstmindis_point
                firstmindis = dis
                firstmindis_point = feap
            else:
                if dis < secondmindis:
                    secondmindis = dis
                    secondmindis_point = feap


        


        
        right_tra, left_tra = [], []
        for i in tra_all:
            if i[2] != firstmindis_point[2] and i[2] != secondmindis_point[2] and i[3] != firstmindis_point[3] and i[3] != \
                    secondmindis_point[3]:
                ratio = cal_ang((i[2], i[3]), (firstmindis_point[2], firstmindis_point[3]),
                                (secondmindis_point[2], secondmindis_point[3]))
                if 0 < ratio < 90:
                    right_tra.append(i)
                if 90 < ratio < 180:
                    left_tra.append(i)

        
        def GetRightTra_Vertex(tra, pointA_x, pointA_y, pointB_x, pointB_y):
            
            maxdis_hori = 0
            maxdis_ver = 0
            slope_AB = (pointB_y - pointA_y) / (pointB_x - pointA_x)
            slope_AC = -(1 / slope_AB)
            b_AC = pointA_y - slope_AC * pointA_x

            for right_tra_p in tra:
                
                p_to_AB = point_to_line_Distance(right_tra_p, (0, 0, pointA_x, pointA_y, 0),
                                                 (0, 0, pointB_x, pointB_y, 0))
                
                p_to_AC = math.fabs(
                    (slope_AC * right_tra_p[2] - right_tra_p[3] + b_AC) / (math.sqrt(1 + math.pow(slope_AC, 2))))

                
                if p_to_AB > maxdis_hori:
                    maxdis_hori = p_to_AB
                if p_to_AC > maxdis_ver:
                    maxdis_ver = p_to_AC

            if maxdis_ver > maxdis_hori:
                square_len = maxdis_ver + 0.001
            else:
                square_len = maxdis_hori + 0.001

            
            delt_a = 1 + math.pow(slope_AC, 2)

            delt_b = -(2 * pointA_x - 2 * slope_AC * (b_AC - pointA_y))

            delt_c = math.pow(b_AC - pointA_y, 2) - math.pow(square_len, 2) + math.pow(pointA_x, 2)

            delt = math.sqrt(delt_b * delt_b - (4 * delt_a * delt_c))

           
            C1_x = (-delt_b + delt) / (2 * delt_a)
            C1_y = slope_AC * C1_x + b_AC

            C2_x = (-delt_b - math.sqrt(delt_b * delt_b - 4 * delt_a * delt_c)) / (2 * delt_a)
            C2_y = slope_AC * C2_x + b_AC

            

            def Get_parallelC_pointRight(slope_AC, slope_AB, C_X, C_Y, dissquare, firstmindis_point_x,
                                         firstmindis_point_y, secondmindis_point_x, secondmindis_point_y):
            
                t1_C1 = slope_AC - slope_AB
                t2_C1 = (slope_AB - slope_AC) * C_X
                t3_C1 = dissquare * dissquare * (1 + slope_AC * slope_AC)
                x1_C2 = ((-(2 * t1_C1 * t2_C1)) + math.sqrt(
                    4 * t1_C1 * t1_C1 * t2_C1 * t2_C1 - 4 * t1_C1 * t1_C1 * (t2_C1 * t2_C1 - t3_C1))) / (
                                    2 * t1_C1 * t1_C1)
                y1_C2 = slope_AB * x1_C2 + C_Y - slope_AB * C_X

                x2_C2 = ((-(2 * t1_C1 * t2_C1)) - math.sqrt(
                    4 * t1_C1 * t1_C1 * t2_C1 * t2_C1 - 4 * t1_C1 * t1_C1 * (t2_C1 * t2_C1 - t3_C1))) / (
                                    2 * t1_C1 * t1_C1)
                y2_C2 = slope_AB * x2_C2 + C_Y - slope_AB * C_X

                ratio_x1_C2 = cal_ang((x1_C2, y1_C2), (firstmindis_point_x, firstmindis_point_y),
                                      (secondmindis_point_x, secondmindis_point_y))
                ratio_x2_C2 = cal_ang((x2_C2, y2_C2), (firstmindis_point_x, firstmindis_point_y),
                                      (secondmindis_point_x, secondmindis_point_y))
                if ratio_x1_C2 < 90:
                    return x1_C2, y1_C2
                if ratio_x2_C2 < 90:
                    return x2_C2, y2_C2

            D1_x, D1_y = Get_parallelC_pointRight(slope_AC, slope_AB, C1_x, C1_y, square_len, pointA_x, pointA_y,
                                                  pointB_x, pointB_y)
            D2_x, D2_y = Get_parallelC_pointRight(slope_AC, slope_AB, C2_x, C2_y, square_len, pointA_x, pointA_y,
                                                  pointB_x, pointB_y)

            return [C1_x, C1_y], [D1_x, D1_y], [D2_x, D2_y], [C2_x, C2_y]

        def GetLeftTra_Vertex(tra, pointA_x, pointA_y, pointB_x, pointB_y):
           
            maxdis_hori = 0
            maxdis_ver = 0
            slope_AB = (pointB_y - pointA_y) / (pointB_x - pointA_x)
            slope_AC = -(1 / slope_AB)
            b_AC = pointA_y - slope_AC * pointA_x

            for right_tra_p in tra:
               
                p_to_AB = point_to_line_Distance(right_tra_p, (0, 0, pointA_x, pointA_y, 0),
                                                 (0, 0, pointB_x, pointB_y, 0))
                
                p_to_AC = math.fabs(
                    (slope_AC * right_tra_p[2] - right_tra_p[3] + b_AC) / (math.sqrt(1 + math.pow(slope_AC, 2))))

               
                if p_to_AB > maxdis_hori:
                    maxdis_hori = p_to_AB
                if p_to_AC > maxdis_ver:
                    maxdis_ver = p_to_AC

            if maxdis_ver > maxdis_hori:
                square_len = maxdis_ver
            else:
                square_len = maxdis_hori

            
            delt_a = 1 + math.pow(slope_AC, 2)

            delt_b = -(2 * pointA_x - 2 * slope_AC * (b_AC - pointA_y))

            delt_c = math.pow(b_AC - pointA_y, 2) - math.pow(square_len, 2) + math.pow(pointA_x, 2)

            delt = math.sqrt(delt_b * delt_b - (4 * delt_a * delt_c))

           
            C1_x = (-delt_b + delt) / (2 * delt_a)
            C1_y = slope_AC * C1_x + b_AC

            C2_x = (-delt_b - math.sqrt(delt_b * delt_b - 4 * delt_a * delt_c)) / (2 * delt_a)
            C2_y = slope_AC * C2_x + b_AC

            

            def Get_parallelC_pointLeft(slope_AC, slope_AB, C_X, C_Y, dissquare, firstmindis_point_x,
                                        firstmindis_point_y, secondmindis_point_x, secondmindis_point_y):
                
                t1_C1 = slope_AC - slope_AB
                t2_C1 = (slope_AB - slope_AC) * C_X
                t3_C1 = dissquare * dissquare * (1 + slope_AC * slope_AC)
                x1_C2 = ((-(2 * t1_C1 * t2_C1)) + math.sqrt(
                    4 * t1_C1 * t1_C1 * t2_C1 * t2_C1 - 4 * t1_C1 * t1_C1 * (t2_C1 * t2_C1 - t3_C1))) / (
                                    2 * t1_C1 * t1_C1)
                y1_C2 = slope_AB * x1_C2 + C_Y - slope_AB * C_X

                x2_C2 = ((-(2 * t1_C1 * t2_C1)) - math.sqrt(
                    4 * t1_C1 * t1_C1 * t2_C1 * t2_C1 - 4 * t1_C1 * t1_C1 * (t2_C1 * t2_C1 - t3_C1))) / (
                                    2 * t1_C1 * t1_C1)
                y2_C2 = slope_AB * x2_C2 + C_Y - slope_AB * C_X

                ratio_x1_C2 = cal_ang((x1_C2, y1_C2), (firstmindis_point_x, firstmindis_point_y),
                                      (secondmindis_point_x, secondmindis_point_y))
                ratio_x2_C2 = cal_ang((x2_C2, y2_C2), (firstmindis_point_x, firstmindis_point_y),
                                      (secondmindis_point_x, secondmindis_point_y))
                if ratio_x1_C2 > 90:
                    return x1_C2, y1_C2
                if ratio_x2_C2 > 90:
                    return x2_C2, y2_C2

            D1_x, D1_y = Get_parallelC_pointLeft(slope_AC, slope_AB, C1_x, C1_y, square_len, pointA_x, pointA_y,
                                                 pointB_x, pointB_y)
            D2_x, D2_y = Get_parallelC_pointLeft(slope_AC, slope_AB, C2_x, C2_y, square_len, pointA_x, pointA_y,
                                                 pointB_x, pointB_y)
            return [C1_x, C1_y], [D1_x, D1_y], [D2_x, D2_y], [C2_x, C2_y]

        rightvec1, rightvec2, rightvec3, rightvec4 = GetRightTra_Vertex(right_tra, firstmindis_point[2],
                                                                        firstmindis_point[3], secondmindis_point[2],
                                                                        secondmindis_point[3])
        leftvec1, leftvec2, leftvec3, leftvec4 = GetLeftTra_Vertex(left_tra, firstmindis_point[2], firstmindis_point[3],
                                                                   secondmindis_point[2], secondmindis_point[3])

        vec_right = [rightvec1, rightvec2, rightvec3, rightvec4]
        vec_left = [leftvec1, leftvec2, leftvec3, leftvec4]

        
        zone_fp=[]

        
        def DivideRecRobust(veccor, point, flag, threshold, robuest_zone):
            
            inAreaPoint = []
            for p in point:
                if isInRectangle(veccor[0][0], veccor[0][1], veccor[1][0], veccor[1][1], veccor[2][0], veccor[2][1],
                                 veccor[3][0], veccor[3][1], p[2], p[3]):
                    inAreaPoint.append(p)

            
            featurepoint_count = 0
            featurepointlist = []
            for feature_point in inAreaPoint:
                # print(feature_point)
                if feature_point[4] == "T":
                    featurepoint_count += 1
                    featurepointlist.append(feature_point)

            
            if featurepoint_count >= threshold:
                
                if flag == True:
                    newveccor1 = [veccor[0], veccor[1],
                                  [(veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2],
                                  [(veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2]]
                    newveccor2 = [[(veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2],
                                  [(veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2], veccor[2],
                                  veccor[3]]
                    newveccor1count, newveccor2count = 0, 0
                    for p in featurepointlist:
                        if isInRectangle(veccor[0][0], veccor[0][1], veccor[1][0], veccor[1][1],
                                         (veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2,
                                         (veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2, p[2],
                                         p[3]):
                            newveccor1count += 1
                        if isInRectangle((veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2,
                                         (veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2,
                                         veccor[2][0], veccor[2][1], veccor[3][0], veccor[3][1], p[2], p[3]):
                            newveccor2count += 1
                    
                    if newveccor1count < threshold and newveccor2count < threshold:
                        robuest_zone.append([inAreaPoint, veccor])
                    else:
                        DivideRecRobust(newveccor1, inAreaPoint, False, threshold, robuest_zone)
                        DivideRecRobust(newveccor2, inAreaPoint, False, threshold, robuest_zone)

                if flag == False:
                    newveccor1 = [veccor[0], [(veccor[0][0] + veccor[1][0]) / 2, (veccor[0][1] + veccor[1][1]) / 2],
                                  [(veccor[2][0] + veccor[3][0]) / 2, (veccor[2][1] + veccor[3][1]) / 2], veccor[3]]
                    newveccor2 = [[(veccor[0][0] + veccor[1][0]) / 2, (veccor[0][1] + veccor[1][1]) / 2], veccor[1],
                                  veccor[2], [(veccor[2][0] + veccor[3][0]) / 2, (veccor[2][1] + veccor[3][1]) / 2]]
                    newveccor1count, newveccor2count = 0, 0
                    for p in featurepointlist:
                        if isInRectangle(veccor[0][0], veccor[0][1], (veccor[0][0] + veccor[1][0]) / 2,
                                         (veccor[0][1] + veccor[1][1]) / 2, (veccor[2][0] + veccor[3][0]) / 2,
                                         (veccor[2][1] + veccor[3][1]) / 2, veccor[3][0], veccor[3][1], p[2], p[3]):
                            newveccor1count += 1
                        if isInRectangle((veccor[0][0] + veccor[1][0]) / 2, (veccor[0][1] + veccor[1][1]) / 2,
                                         veccor[1][0], veccor[1][1], veccor[2][0], veccor[2][1],
                                         (veccor[2][0] + veccor[3][0]) / 2, (veccor[2][1] + veccor[3][1]) / 2, p[2],
                                         p[3]):
                            newveccor2count += 1
                    
                    if newveccor1count < threshold and newveccor2count < threshold:
                        robuest_zone.append([inAreaPoint, veccor])
                    else:
                        DivideRecRobust(newveccor1, inAreaPoint, True, threshold, robuest_zone)
                        DivideRecRobust(newveccor2, inAreaPoint, True, threshold, robuest_zone)
            else:
                if len(inAreaPoint) > 0:
                    robuest_zone.append([inAreaPoint, veccor])

        def area_of_quadrangle(x1, y1, x2, y2, x3, y3, x4, y4):
            
            area1 = abs(0.5 * (x1 * y2 + x2 * y3 + x3 * y4 + x4 * y1 - x2 * y1 - x3 * y2 - x4 * y3 - x1 * y4))
            area2 = abs(0.5 * (x1 * y3 + x3 * y4 + x4 * y2 + x2 * y1 - x3 * y1 - x4 * y2 - x2 * y3 - x1 * y4))
            
            area = area1 + area2
            return area

        def DivedeRecFragile(veccor, point, flag, threshold, zone):
            inAreaPoint = []
            for p in point:
                if isInRectangle(veccor[0][0], veccor[0][1], veccor[1][0], veccor[1][1], veccor[2][0], veccor[2][1],
                                 veccor[3][0], veccor[3][1], p[2], p[3]):
                    inAreaPoint.append(p)

            
            recarea = area_of_quadrangle(veccor[0][0], veccor[0][1], veccor[1][0], veccor[1][1], veccor[2][0],
                                         veccor[2][1], veccor[3][0], veccor[3][1])
            # print(recarea)

            
            if recarea >= threshold and len(inAreaPoint) > 0:
                
                if flag == True:
                    newveccor1 = [veccor[0], veccor[1],
                                  [(veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2],
                                  [(veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2]]
                    newveccor2 = [[(veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2],
                                  [(veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2], veccor[2],
                                  veccor[3]]

                    a1 = area_of_quadrangle(veccor[0][0], veccor[0][1], veccor[1][0], veccor[1][1],
                                            (veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2,
                                            (veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2)

                    a2 = area_of_quadrangle((veccor[0][0] + veccor[3][0]) / 2, (veccor[0][1] + veccor[3][1]) / 2,
                                            (veccor[1][0] + veccor[2][0]) / 2, (veccor[1][1] + veccor[2][1]) / 2,
                                            veccor[2][0], veccor[2][1], veccor[3][0], veccor[3][1])

                    
                    if a1 < threshold and a2 < threshold:
                        zone.append(inAreaPoint)
                    else:
                        DivedeRecFragile(newveccor1, inAreaPoint, False, threshold, zone)
                        DivedeRecFragile(newveccor2, inAreaPoint, False, threshold, zone)

                if flag == False:
                    newveccor1 = [veccor[0], [(veccor[0][0] + veccor[1][0]) / 2, (veccor[0][1] + veccor[1][1]) / 2],
                                  [(veccor[2][0] + veccor[3][0]) / 2, (veccor[2][1] + veccor[3][1]) / 2], veccor[3]]
                    newveccor2 = [[(veccor[0][0] + veccor[1][0]) / 2, (veccor[0][1] + veccor[1][1]) / 2], veccor[1],
                                  veccor[2], [(veccor[2][0] + veccor[3][0]) / 2, (veccor[2][1] + veccor[3][1]) / 2]]

                    a3 = area_of_quadrangle(veccor[0][0], veccor[0][1], (veccor[0][0] + veccor[1][0]) / 2,
                                            (veccor[0][1] + veccor[1][1]) / 2, (veccor[2][0] + veccor[3][0]) / 2,
                                            (veccor[2][1] + veccor[3][1]) / 2, veccor[3][0], veccor[3][1])

                    a4 = area_of_quadrangle((veccor[0][0] + veccor[1][0]) / 2, (veccor[0][1] + veccor[1][1]) / 2,
                                            veccor[1][0], veccor[1][1], veccor[2][0], veccor[2][1],
                                            (veccor[2][0] + veccor[3][0]) / 2, (veccor[2][1] + veccor[3][1]) / 2)

                    
                    if a3 < threshold and a4 < threshold:
                        zone.append(inAreaPoint)
                    else:
                        DivedeRecFragile(newveccor1, inAreaPoint, True, threshold, zone)
                        DivedeRecFragile(newveccor2, inAreaPoint, True, threshold, zone)
            else:
                if len(inAreaPoint) > 0:
                    zone.append(inAreaPoint)

        
        DivideRecRobust(vec_right, right_tra, True, robust_threshold,zone_fp)
        DivideRecRobust(vec_left, left_tra, True, robust_threshold,zone_fp)

        final_zone = []
        areathreshold = 0.00001

        for z in zone_fp:
            zone_tmp = z[0]

            count = 0
            for i in zone_tmp:
                if i[4] == "T":
                    count = count + 1
            # print(count)
            if robust_threshold < count:
                DivedeRecFragile(z[1], zone_tmp, True, fragile_threshold, final_zone)
            else:
                final_zone.append(zone_tmp)


        

        
        embed_bit = str2bitarray(embedbit).tolist()

        slope_AB = (secondmindis_point[3] - firstmindis_point[3]) / (secondmindis_point[2] - firstmindis_point[2])
        b_AB = firstmindis_point[3] - (slope_AB * firstmindis_point[2])

        

        for zone in final_zone:
            changedata = copy.deepcopy(zone)

            
            Selection_sort_bydis(changedata, firstmindis_point, secondmindis_point)


            fpcount = 0
            for fp in zone:
                if fp[4] == "T":
                    fpcount += 1



            
            if fpcount >= robust_threshold:
                ncvalue=RoubustWaterMarkEmbed(changedata, embed_bit, 2, 1,firstmindis_point,secondmindis_point)

            
            fragile_result=FragileWaterMarkEmbed(changedata, "ID", 2, 1,firstmindis_point,secondmindis_point)



