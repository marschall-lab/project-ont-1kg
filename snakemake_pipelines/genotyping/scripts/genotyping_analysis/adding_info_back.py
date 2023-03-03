import sys
import cyvcf2

def read_input_header(vcf):
    info = []
    for i in vcf.header_iter():
        if i["HeaderType"] == "INFO":
             info.append(i)
    return info

def update_vcf(vcf, info):
    for i in info:
        vcf.add_info_to_header(i)
    return vcf

callset_vcf = sys.argv[1] #This is the VCF that will be changed
vcf = sys.argv[2] #This is the VCF whose INFO field will be added to callset_vcf
out = callset_vcf[0:-7]+"_info_added.vcf"

info_vcf = cyvcf2.VCF(vcf)
info_iter = iter(info_vcf)
callset_vcf = cyvcf2.VCF(callset_vcf)
callset_iter = iter(callset_vcf)
header = read_input_header(info_vcf)
info_dict = []
for h in header:
    d = {"Type": h["Type"], "Number": h["Number"], "ID": h["ID"], "Description": h["Description"][1:-1]}
    info_dict.append(d)

vcf_list = []
callset_vcf = update_vcf(callset_vcf, info_dict)
vcf_out = cyvcf2.Writer(out, callset_vcf)

done_looping = False

info_list = ["CONFLICT", "LV", "PS", "AT", "ID"]
while not done_looping:
    try:
        var_info = next(info_iter)
        info = []
        for i in info_list:
            try:
                info.append(var_info.INFO[i])
            except KeyError:
                info.append(-1)
        var_callset = next(callset_iter)
        assert (var_callset.POS == var_info.POS)
        for i in range(len(info)):
            if info[i] == -1:
                continue
            try:
                var_callset.INFO[info_list[i]] = info[i]
            except AttributeError:
                var_callset.INFO[info_list[i]] = ",".join([str(x) for x in info[i]])
        vcf_out.write_record(var_callset)
        
    except StopIteration:
        done_looping = True
        
callset_vcf.close()
info_vcf.close()
vcf_out.close()