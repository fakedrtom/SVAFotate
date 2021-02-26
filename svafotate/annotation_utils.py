import os
import sys
from .utils import get_feature



def create_max_annotations(my_list,vcf,req_sources):
    ## for each item in the list
    ## create Max_AF,Max_Het,Max_HomAlt annotations
    for i in my_list:
        af = ["Max",str(i),"AF"]
        vcf.add_info_to_header({"ID": "_".join(af), "Description": "The maximum " + str(i) + " AF from all matching SVs across all specified data sources (" + ",".join(req_sources) + ")", "Type": "Float", "Number": "1"})
        het = ["Max",str(i),"Het"]
        vcf.add_info_to_header({"ID": "_".join(het), "Description": "The maximum " + str(i) + " Het count from all matching SVs across all specified data sources (" + ",".join(req_sources) + ")", "Type": "Integer", "Number": "1"})
        homalt = ["Max",str(i),"HomAlt"]
        vcf.add_info_to_header({"ID": "_".join(homalt), "Description": "The maximum " + str(i) + " HomAlt count from all matching SVs across all specified data sources (" + ",".join(req_sources) + ")", "Type": "Integer", "Number": "1"})


def create_best_annotations(source,my_list,vcf):
    ## for each item in the list
    ## create Best_source_AF,Best_source_Het,Best_source_HomAlt annotations
    for i in my_list:
        af = ["Best",str(source),str(i),"AF"]
        vcf.add_info_to_header({"ID": "_".join(af), "Description": "The best " + str(i) + " AF match for " + str(source), "Type": "Float", "Number": "1"})
        het = ["Best",str(source),str(i),"Het"]
        vcf.add_info_to_header({"ID": "_".join(het), "Description": "The best " + str(i) + " Het count match for " + str(source), "Type": "Integer", "Number": "1"})
        homalt = ["Best",str(source),str(i),"HomAlt"]
        vcf.add_info_to_header({"ID": "_".join(homalt), "Description": "The best " + str(i) + " HomAlt count match " + str(source), "Type": "Integer", "Number": "1"})


#def join_annotations(source):
    ## create annotations specific to sources requested
    ## makes "best" annotations corresponding to specified fields
#    vcf.add_info_to_header({"ID": "Best_" + source, "Description": "The " + source + " ID of the best matching SV", "Type": "String", "Number": "1"})
#    for i in source_cols[source]:
#        annotation = bed_headers[i]
#        annotation_type = header_types[annotation]
#        vcf.add_info_to_header({"ID": "Best_" + source + "_" + annotation, "Description": "The best match " + annotation + " for " + source, "Type": annotation_type, "Number": "1"})


#def add_best_annotations(source,my_list):
    ## writes out best annotations to vcf
#    if len(my_list) > 1:
#        print(source + " has too many best")
#    else:
#        v.INFO["Best_" + source] = my_list[0]
#        for i in source_cols[source]:
#            annotation = bed_headers[i]
#            annotation_type = header_types[annotation]
#            v.INFO["Best_" + source + "_" + annotation] = datas[source][my_list[0]][i]


def write_max_values(sv_id,my_dict,features,req_sources,datas,header_cols,header_types,v):
    ## features is a list of annotations
    ## finds max values for list of annotations
    ## adds those max values to output VCF
    for i in features:
        anno = "Max_" + str(i)
        max_val = float(0)
        all_vals = []
        for source in req_sources:
            matches_vals = []
            if sv_id in my_dict[source]:
                matches_vals = get_feature(source,my_dict[source][sv_id],header_cols[i],datas)
            all_vals.extend(matches_vals)
        if len(all_vals) > 0:
            if header_types[i] == "Float":
                max_val = max(list(map(float,all_vals)))
            if header_types[i] == "Integer":
                max_val = max(list(map(int,all_vals)))
        v.INFO[anno] = max_val


def write_best_values(source,sv_id,my_dict,features,datas,header_cols,v):
    ## features is a list of annotations
    ## finds best values for list of annotations
    ## adds those best values to output VCF
    for i in features:
        anno = "Best_" + str(source) + "_" + str(i)
        best_val = []
        if sv_id in my_dict[source]:
            matches_vals = get_feature(source,my_dict[source][sv_id],header_cols[i],datas)
            best_val.extend(matches_vals)
        if len(best_val) == 0:
            best_val = [0]
        v.INFO[anno] = best_val[0]


