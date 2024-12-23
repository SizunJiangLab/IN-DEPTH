#################################################################################
#
# Merge ROI
#
#################################################################################

import pandas as pd

f = open(r"merged_24_25_raw.csv")
colname = f.readline()
cell_list = ["CD4Treg","CD8T","CD4T","M1","Myeloid","BCL6pB","BCL6nB","M2"]
roi_list = []
roi_dict = {}
sampler_list = colname.strip().split(",")
flag = 0

for item in sampler_list:
    if len(item) > 3:
        roi_name = item.strip().split(".")[3]
        celltype = item.strip().split(".")[6].strip('/"')
        #roi_name = item.strip().split("|")[1].strip(" ")
        if celltype in cell_list:
            print(celltype)
            if roi_name not in roi_dict.keys():
                roi_dict[roi_name] = [flag]
            else:
                roi_dict[roi_name].append(flag)
            if roi_name not in roi_list:
                roi_list.append(roi_name)
        flag += 1
df = pd.read_csv(r"merged_24_25_raw.csv",index_col=0)

for num in roi_list:
    tmp_lst = roi_dict[num]
    df[num] = df.iloc[:,tmp_lst].sum(axis=1)
#    df[num] = df.iloc[:,tmp_lst].mean(axis=1)
df = df.iloc[:,332:]
df.to_csv("roi_merge_data_raw_clean.csv")
f.close()

#################################################################################
#
# Calculate pearson correlation and coefficient, RMS and JSD
#
#################################################################################

# ground_truth.txt sc-RNA
# CD4T CD4 T(9)
# Myeloid myeloid(0)
# M1 M1 Macrophages(3)
# CD8T CD8 T(7)
# DC DC(2) 
# M2 M2 Macrophages(1)
# BCL6pB GCBC
# BCL6nB MBC
import math
f_out = open(r"indepth_result.txt","w")
f_out2 = open(r"indepth_result2.txt","w")
f_out3 = open(r"indepth_result3.txt","w")
match_dict = {"CD4T": 5, "Myeloid": 3, "M1": 1, "CD8T": 4, "DC": 2, "M2": 0, "CD4Treg": 8, "BCL6nB": 6, "BCL6pB": 7} #dtangle
#match_dict = {"CD4T": 5, "Myeloid": 3, "CD8T": 4, "DC": 2, "M2": 0, "BCL6nB": 6, "BCL6pB": 7, "CD4Treg": 8, "M1": 1} #MuSiC
#match_dict = {"CD4T": 1, "Myeloid": 16, "M1": 6, "CD8T": 2, "DC": 3, "M2": 7} #stereoscope
#match_dict = {"CD4T": 4, "Myeloid": 2, "CD8T": 3, "DC": 1, "M2": 0, "BCL6nB": 5, "BCL6pB": 6, "CD4Treg": 7, "M1": 8} #Cybersort
#match_dict = {"CD4T": 4, "Myeloid": 2, "CD8T": 3, "DC": 1, "M2": 0, "BCL6nB": 5, "BCL6pB": 6, "CD4Treg": 7, "M1": 8} #spatial_decon
match_dict = {"CD4T": 0, "Myeloid": 8, "CD8T": 1, "DC": 2, "M2": 5, "BCL6nB": 6, "BCL6pB": 3, "CD4Treg": 7, "M1": 4} #stereoscope


predict_dict = {}
f = open(r"dtangle_result.tsv")
for line in f:
    if line[1] == "X" or line[0] == "X":
        tmp_list = line.strip().split("\t")
        #tmp_list = line.strip().split(',')
        predict_dict[tmp_list[0].strip('/"').strip("X")] = []
        for i in range(1,len(tmp_list),1):
            predict_dict[tmp_list[0].strip('/"').strip("X")].append(float(tmp_list[i]))
f.close() 

ground_truth_dict = {}
f2 = open(r"ground_truth.txt")
for line in f2:
    if line.startswith("ROI"):
        roi_name = line.strip().split(" ")[1]
        ground_truth_dict[roi_name] = {}
    else:
        cell_type_name = line.strip().split(" ")[0].strip(":")
        true_ratio = float(line.strip().split(" ")[1])
        ground_truth_dict[roi_name][cell_type_name] = true_ratio
f2.close()

sum_ratio = []
flag = 0
for key in ground_truth_dict.keys():
    x = []
    y = []
    E_xy = 0
    E_x = 0
    E_y = 0
    sigma_x = 0
    sigma_y = 0
    for key2 in match_dict.keys():
        # if key2 in ground_truth_dict[key].keys():
        #     if key2 == "CD4T":
        #         x.append(ground_truth_dict[key][key2])
        #         tmp_sum = 0
        #         for j in match_dict[key2]:
        #             tmp_sum += predict_dict[key][j]
        #         y.append(tmp_sum)
        #     else:
        #         x.append(ground_truth_dict[key][key2])
        #         y.append(predict_dict[key][match_dict[key2]])
        if key2 in ground_truth_dict[key].keys():
                x.append(ground_truth_dict[key][key2])
                y.append(predict_dict[key][match_dict[key2]])
    for i in range(len(x)):
        E_xy += x[i]*y[i]
        E_x += x[i]
        E_y += y[i]
    if len(x) == 0:
        E_xy = 0 
        E_x = 0
        E_y = 0
        print(f"When {key} len(x) == 0")
    else:
        E_xy /= len(x)
        E_x /= len(x)
        E_y /= len(x)
        flag += 1
    for i in range(len(x)):
        sigma_x += (x[i] - E_x)*(x[i] - E_x)
        sigma_y += (y[i] - E_y)*(y[i] - E_y)
    if len(x) == 0:
        sigma_x = 0
        sigma_y = 0
    else:
        sigma_x /= len(x)
        sigma_y /= len(x)
    sigma_x = sigma_x**0.5
    sigma_y = sigma_y**0.5
    if (sigma_x == 0) or (sigma_y == 0):
        #sum_ratio += 0
        sum_ratio.append(0)
    else:  
        #sum_ratio += (E_xy - E_y*E_x)/(sigma_x*sigma_y)
        sum_ratio.append((E_xy - E_y*E_x)/(sigma_x*sigma_y))
#sum_ratio /= flag
#print(f"The Pearson correlation coefficient is:{sum_ratio}")
flag = 0
for key in ground_truth_dict.keys():
    #f_out.write(f"ROI{key}: {sum_ratio[flag]}\n")
    f_out.write(f"{sum_ratio[flag]}\n")
    flag += 1

sum_ratio = []
flag = 0
for key in ground_truth_dict.keys():
    x = []
    y = []
    diff = 0
    for key2 in match_dict.keys():
        if key2 in ground_truth_dict[key].keys():
            x.append(ground_truth_dict[key][key2])
            y.append(predict_dict[key][match_dict[key2]])
        # if key2 in ground_truth_dict[key].keys():
        #     if key2 == "CD4T":
        #         x.append(ground_truth_dict[key][key2])
        #         tmp_sum = 0
        #         for j in match_dict[key2]:
        #             tmp_sum += predict_dict[key][j]
        #         y.append(tmp_sum)
        #     else:
        #         x.append(ground_truth_dict[key][key2])
        #         y.append(predict_dict[key][match_dict[key2]])
    for i in range(len(x)):
        diff += (x[i] - y[i])**2
    if len(x) == 0:
        diff = 0 
        print(f"When {key} len(x) == 0")
    else:
        diff /= len(x)
        flag += 1
    diff = diff**0.5
    #sum_ratio += diff
    sum_ratio.append(diff)
#sum_ratio /= flag
#print(f"The RMS is:{sum_ratio}")
flag = 0
for key in ground_truth_dict.keys():
    #f_out.write(f"ROI{key}: {sum_ratio[flag]}\n")
    f_out2.write(f"{sum_ratio[flag]}\n")
    flag += 1

sum_ratio3 = []
flag = 0
for key in ground_truth_dict.keys():
    x = []
    y = []
    d1 = 0
    d2 = 0
    tmp_sum1 = 0
    tmp_sum2 = 0
    for key2 in match_dict.keys():
        if key2 in ground_truth_dict[key].keys():
            x.append(ground_truth_dict[key][key2])
            y.append(predict_dict[key][match_dict[key2]])
        # if key2 in ground_truth_dict[key].keys():
        #     if key2 == "CD4T":
        #         x.append(ground_truth_dict[key][key2])
        #         tmp_sum = 0
        #         for j in match_dict[key2]:
        #             tmp_sum += predict_dict[key][j]
        #         y.append(tmp_sum)
        #     else:
        #         x.append(ground_truth_dict[key][key2])
        #         y.append(predict_dict[key][match_dict[key2]])
    for i in range(len(x)):
        tmp_sum1 += x[i]
        tmp_sum2 += y[i]
    x.append(1 - tmp_sum1)
    y.append(1 - tmp_sum2)
    if ((1 - tmp_sum1) < 0) and ((1 - tmp_sum2) < 0):
        raise Exception(f"Some thing wrong in ROI: {key}")
    for i in range(len(x)):
        if x[i] == 0 or y[i] == 0 or x[i] <= 0 or y[i] <= 0:
            d1 += 0
            d2 += 0
        else:
            d1 += x[i]*math.log2(x[i]/(x[i]+y[i]))
            d2 += y[i]*math.log2(y[i]/(x[i]+y[i]))
    flag += 1
    sum_ratio3.append(0.5*d1 + 0.5 *d2 + 1)
    #sum_ratio3 += 0.5*d1 + 0.5 *d2 + 1
#sum_ratio3 /= flag
#print(f"The Jensen-Shannon divergence is:{sum_ratio3}")
#f_out.write(f"The Jensen-Shannon divergence is:\n")
flag = 0
for key in ground_truth_dict.keys():
    #f_out.write(f"ROI{key}: {sum_ratio[flag]}\n")
    f_out3.write(f"{sum_ratio3[flag]}\n")
    flag += 1

f_out.close()
f_out2.close()
f_out3.close()




















