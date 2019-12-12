__author__ = 'carlos.company@crg.eu'
# -*- coding utf-8 -*-

#MODULES
import sys
import re
import optparse

#BODY FUNTIONS
def options_arg():
    usage = "usage: %prog -p <TAX LEVEL> -t <TAXONOMY FILE e.g. desc.an.cons.taxonomy> -s <ABUNDANCE FILE e.g. desc.an.0.03.subsample.relabund > [OP] -w <OUTPUT FILE>"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-p', '--tax-level', help='Silva Taxonomy level you want to compare', dest="level" )
    parser.add_option('-t', '--taxonomy', help='Taxonomy File (OTUs-taxonomy)', dest="tax" )
    parser.add_option('-s', '--samples', help='Tax File', dest="samples" )
    parser.add_option('-w', '--write', action='store_true',help='ouput text File', dest="wotus" )
    (opts,args) = parser.parse_args()
    if opts.level and opts.tax and opts.samples:pass
    else: parser.print_help()
    return (opts)
def __main__ ():
    data_tax_requested = tax_level()
    tax_sum_data = otus_count(data_tax_requested)
    end_table(tax_sum_data)

#AUXILIAR MODULES
def write_file(data2write):
    name_file = 'mothur_count_taxlevel'+str(tax_num)+'.txt'
    f = open(name_file,'w')
    for i in data2write:
        a = i + '\n'
        f.write(a)
    f.close()
def tax_level():
    line_file = 0
    list_data = {}
    list_data_f = {}
    for i in file_tax:
        if line_file > 0:
            i = i.strip('\n'); i= i.split('\t')
            tax_id = i[2].split(';')
            tax_filt = tax_id[0:tax_num]; tax_filt=';'.join(tax_filt)
            tax_filt = tax_filt.translate(None,'\"\(\)0123456789')
            if((len(tax_id)-1) == tax_num ):
                list_data[i[0]] =  tax_filt
                list_data_f[i[0]] = [float(i[1]), tax_filt]
        line_file +=1
    end_del = [list_data, list_data_f]
    return(end_del)

def samples_dict():
    #samples_data = open('desc.an.0.03.subsample.relabund')
    line_f = 0
    names_otus = []
    values_otus = {}
    for i in samples_data:
        i = i.strip('\n')
        aux_i = i.split('\t')
        if line_f <= 0:
            names_otus = aux_i[3:len(aux_i)-1]
        elif line_f > 0:
            values_otus[aux_i[1]] = aux_i[3:len(aux_i)-1]
        else: pass
        line_f +=1
    res_table = {}
    line_f = 0
    for i in names_otus:
        indv_data = {}
        line_sub = 0
        for p in values_otus.keys():
            #print i + '\t' + str(p) + '\t' + str(values_otus[p][line_f])
            indv_data[p] = str(values_otus[p][line_f])
            line_sub +=1
        res_table[i] = indv_data
        line_f +=1
    return (res_table)
def otus_count(dict_data):
    data_tax = set(dict_data[0].values())
    dict_key = dict_data[0].keys()
    pos_data_tax = {}
    del_res = {}
    for tax in data_tax:
        pos_data_tax[tax] = []
        for i in dict_key:
            if dict_data[0][i] == tax: pos_data_tax[tax].append(i)

    #Sum by Sample
    #for tax in pos_data_tax.keys():
    list_elements = {}
    data_sample = samples_dict()
    for tax in pos_data_tax.keys():
        aux_list = pos_data_tax[tax]
        list_elements[tax] = {}
        for i in data_sample[data_sample.keys()[1]].keys():
            list_elements[tax][i] = float(0)

        for i in aux_list:
            if data_sample.has_key(i):
                for k in data_sample[i].keys():
                    list_elements[tax][k] +=  float(data_sample[i][k])
                    #print tax + '\t' + i + '\t' + k + '\t' + str(data_sample[i][k])

    #Sum GLOBAL Taxonomy Parts
    for tax in pos_data_tax.keys():
        aux_list = pos_data_tax[tax]
        sum_total = 0
        for i in aux_list:
            sum_total += dict_data[1][i][0]
        del_res[tax] = float(sum_total)

    end_res = [del_res,list_elements]
    return (end_res)
def end_table(end_list_data):
    header_file =  'Taxon_Level_'+ str(tax_num)+ '\t' + 'Tax_Number' + '\tS' + "\tS".join(map(str,end_list_data[1].values()[0].keys()))+'\t'+'Total'
    print header_file
    data_write = [header_file]
    for i in end_list_data[1].keys():
        a = end_list_data[1][i].values()
        end_id=i.split(';')
        end_id=end_id[tax_num-1]
        b= end_id + '\t' + i + '\t' +"\t".join(map(str, a)) + '\t' + str(end_list_data[0][i])
        print b
        data_write.append(b)
    if opts.wotus:
        write_file(data_write)
    else:pass

#Calling
try:
    opts = options_arg()
    if opts.level and opts.tax and opts.samples:
        tax_num = int(opts.level)
        file_tax = open(opts.tax)
        samples_data = open(opts.samples)
        __main__()
    else:
        parser.print_help()
except:pass

