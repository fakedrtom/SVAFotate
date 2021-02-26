import os
import sys



def process_bed_source(bed_file,covAF,uniqAF,size_limit):
    ## process bed file containing SV AFs
    ## save info into datas in source specific manner
    ## create source specific lists for use in pyranges

    from collections import defaultdict
    import gzip
    import io

    sources = set()
    datas = defaultdict(lambda: defaultdict(list))
    bed_lists = defaultdict(lambda: defaultdict(list))
    bed_headers = []
    cov_lists = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    uniq_lists = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    bed = gzip.open(bed_file, "rt", encoding="utf-8") if bed_file.endswith(".gz") else io.open(bed_file, "rt", encoding="utf-8")
    print("\nReading the following data source file: {}".format(bed_file))

    header = []
    header_found = False
    features = ["chrom","start","end","svlen","svtype","sv_id"]
    for line in bed:

        ## get header
        if line.startswith("#"):
            header = line.strip().replace("#","").split("\t")
            header_found = True
            bed_headers.append(header[4])
            bed_headers.extend(header[7:len(header)])
            continue

        ## Check for header
        if header_found == False:
            raise ValueError("A header is required in the bed source file")

        ## Set up data
        line_dict = dict(zip(header, line.strip().split("\t")))
        source = line_dict["SOURCE"]
        sv_id = line_dict["SV_ID"]
        svlen = line_dict["SVLEN"]
        svtype = line_dict["SVTYPE"]
        start = line_dict["START"]
        end = line_dict["END"]
        size = int(end) - int(start)

        ## Collect info
        datas[source][sv_id].append(line_dict["SVTYPE"])
        datas[source][sv_id].extend([line_dict[field] for field in header[7:len(header)]])
    
        ## define AF for each line
        ## MCNV from gnomAD currently lists lowAF:highAF
        ## correct these AFs to just the highAF
        af = line_dict["AF"]
        if ':' in af:
            af = af.split(':')[1]

        ## add freature info for
        for i,key in enumerate(features):
            bed_lists[source][features[i]].append(line_dict[key.upper()])
            if float(af) > covAF and size <= size_limit:
                cov_lists[source][svtype][features[i]].append(line_dict[key.upper()])
            if float(af) > uniqAF and size <= size_limit:
                uniq_lists[source][svtype][features[i]].append(line_dict[key.upper()])

        ## add sources sources 
        sources.add(source)

    bed.close()

    return(sources,datas,bed_lists,bed_headers,cov_lists,uniq_lists)


def process_pickled_source(pickled_source,coAF,uniqAF,size_limit):

    import pickle

    print("\nReading pickled data source")
    pickle_file = open(pickled_source,"rb")
    pickled_data = pickle.load(pickle_file)
    pickle_file.close()
    
    return(pickled_data["sources"], pickled_data["datas"], pickled_data["bed_lists"], pickled_data["bed_headers"], pickled_data["cov_lists"], pickled_data["uniq_lists"])


def get_feature(source,my_list,col,datas):
    ## takes a list of SV_ID_b and returns the feature for the specified col
    ## returns list of requested feature
    features = []
    for i in my_list:
        feature = datas[source][i][col]
        if feature != "NA": 
            features.append(feature)
    return(features)
