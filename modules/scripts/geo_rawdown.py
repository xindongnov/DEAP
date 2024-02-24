#!/usr/bin/env python

# ================================
# @auther: Xin Dong
# @email: xindong9511@gmail.com
# @date: Sep 2019
# ================================

import os, sys
import urllib.request
import re
import argparse
import subprocess

def gzip_fastq(lay_type,path,gsm):
    global compress
    if compress == True:
        if lay_type == "SINGLE":
            os.system('gzip %s/%s.fastq \n' % (path, gsm))
        elif lay_type == 'PAIRED':
            os.system('gzip %s/%s_R1.fastq \n'%(path, gsm))
            os.system('gzip %s/%s_R2.fastq \n'%(path, gsm))
    else:
        return None

def gunzip_fastq(lay_type,path,gsm):
    global compress
    if compress != True:
        if lay_type == "SINGLE":
            os.system('gunzip %s/%s.fastq.gz \n' % (path, gsm))
        elif lay_type == 'PAIRED':
            os.system('gunzip %s/%s_R1.fastq.gz \n'%(path, gsm))
            os.system('gunzip %s/%s_R2.fastq.gz \n'%(path, gsm))
    else:
        return None

def checkSRA(fsra):
    try:
        f = os.path.basename(fsra)
        cmd_output = subprocess.getoutput('vdb-validate %s' % fsra)
        if ("'%s' is consistent" % f in cmd_output) or ("'%s' metadata: md5 ok" % f in cmd_output):
            sys.stdout.write("sra file OK: %s\n" % fsra)
            return True
        else:
            sys.stdout.write("sra file problem: %s\n" % fsra)
            return False
    except:
        return False

def checkFastqSize(fastq):
    if os.path.getsize(fastq) < 100:
        os.system("rm %s" % fastq)
        return False
    else:
        return True

def downloadFastqByPrefetch(path,gsm,srr,lay_type):
    cat_file1 = ''
    cat_file2 = ''
    if lay_type == 'SINGLE':
        for i in range(len(srr)):
            fsra = '%s/%s_temp%s.sra'%(path, gsm, i+1)
            os.system('prefetch %s -o %s \n'%(srr[i], fsra))
            if os.path.exists(fsra) != True:
                return False
            if checkFastqSize(fsra) != True:
                return False
            if checkSRA(fsra) != True:
                return False
            os.system('echo "+++fastq-dump++++"')
            os.system('\nfastq-dump %s/%s_temp%s.sra -O %s \n' % (path, gsm, i+1, path))
            cat_file1 = cat_file1 + '%s/%s_temp%s.fastq ' % (path, gsm,i+1) 
        os.system('cat %s> %s/%s.fastq \n'%(cat_file1,path,gsm))
        os.system('rm %s/%s_temp* \n'%(path, gsm))
        gzip_fastq(lay_type,path,gsm)
    elif lay_type == 'PAIRED':
        for i in range(len(srr)):
            fsra = '%s/%s_temp%s.sra'%(path, gsm, i+1)
            os.system('prefetch %s -o %s \n'%(srr[i], fsra))
            if os.path.exists(fsra) != True:
                return False
            if checkFastqSize(fsra) != True:
                return False
            os.system('echo "+++fastq-dump++++" \n')
            os.system('\nfastq-dump --split-files %s/%s_temp%s.sra -O %s \n'%(path, gsm,i+1, path))
            cat_file1 = cat_file1 + '%s/%s_temp%s_1.fastq '%(path, gsm,i+1)
            cat_file2 = cat_file2 + '%s/%s_temp%s_2.fastq '%(path, gsm,i+1)
        os.system('cat %s> %s/%s_R1.fastq \n'%(cat_file1, path,gsm))
        os.system('cat %s> %s/%s_R2.fastq \n'%(cat_file2, path,gsm))
        os.system('rm %s/%s_temp* \n'%(path, gsm))
        gzip_fastq(lay_type,path,gsm)
    return True


def downloadFastqFromEBI(path,gsm,srr,lay_type):
    ebi_ftp = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
    if lay_type == "SINGLE":
        cat_file = ""
        if len(srr) == 1:
            s = srr[0]
            down_link = ebi_ftp + s[0:6] + "/00" + s[-1] + "/" + s + "/" + s + ".fastq.gz"
            temp_fq = '%s/%s.fastq.gz'% (path, gsm)
            os.system('wget %s -O %s \n'%(down_link, temp_fq))
            if checkFastqSize(temp_fq) != True:
                return False
            gunzip_fastq(lay_type,path,gsm)
        else:
            for i in range(len(srr)):
                down_link = ebi_ftp + srr[i][0:6] + "/00" + srr[i][-1] + "/" + srr[i] + "/" + srr[i] + ".fastq.gz"
                temp_fq = '%s/%s_temp%s.fastq.gz'% (path, gsm, i+1)
                os.system('wget %s -O %s \n'%(down_link, temp_fq))
                if checkFastqSize(temp_fq) != True:
                    return False
                os.system('gunzip %s \n'%(temp_fq))
                cat_file = cat_file + '%s/%s_temp%s.fastq '%(path, gsm,i+1)
            os.system('cat %s> %s/%s.fastq \n'%(cat_file, path, gsm))
            os.system('rm %s/%s_temp* \n'%(path, gsm))
            gzip_fastq(lay_type,path,gsm)
            gunzip_fastq(lay_type,path,gsm)
        return True
    elif lay_type == "PAIRED":
        cat_file1 = ""
        cat_file2 = ""
        if len(srr) == 1:
            s = srr[0]
            down_link = ebi_ftp + s[0:6] + "/00" + s[-1] + "/" + s + "/" + s + "_1.fastq.gz"
            temp_fq = '%s/%s_R1.fastq.gz'% (path, gsm)
            os.system('wget %s -O %s \n'%(down_link, temp_fq))
            if checkFastqSize(temp_fq) != True:
                return False
            gunzip_fastq(lay_type,path,gsm)
            down_link = ebi_ftp + s[0:6] + "/00" + s[-1] + "/" + s + "/" + s + "_2.fastq.gz"
            temp_fq = '%s/%s_R2.fastq.gz'% (path, gsm)
            os.system('wget %s -O %s \n'%(down_link, temp_fq))
            if checkFastqSize(temp_fq) != True:
                return False
            gunzip_fastq(lay_type,path,gsm)
        else:
            for i in range(len(srr)):
                # Read 1
                down_link = ebi_ftp + srr[i][0:6] + "/00" + srr[i][-1] + "/" + srr[i] + "/" + srr[i] + "_1.fastq.gz"
                temp_fq = '%s/%s_temp%s_1.fastq.gz'% (path, gsm, i+1)
                os.system('wget %s -O %s \n'%(down_link, temp_fq))
                if checkFastqSize(temp_fq) != True:
                    return False
                os.system('gunzip %s \n'%(temp_fq))
                cat_file1 = cat_file1 + '%s/%s_temp%s_1.fastq '%(path, gsm,i+1)
                # Read 2
                down_link = ebi_ftp + srr[i][0:6] + "/00" + srr[i][-1] + "/" + srr[i] + "/" + srr[i] + "_2.fastq.gz"
                temp_fq = '%s/%s_temp%s_2.fastq.gz'% (path, gsm, i+1)
                os.system('wget %s -O %s \n'%(down_link, temp_fq))
                if checkFastqSize(temp_fq) != True:
                    return False
                gzip_fastq(lay_type,path,gsm)
                cat_file2 = cat_file2 + '%s/%s_temp%s_2.fastq '%(path, gsm,i+1)
            os.system('cat %s> %s/%s_R1.fastq \n'%(cat_file1, path, gsm))
            os.system('cat %s> %s/%s_R2.fastq \n'%(cat_file2, path, gsm))
            os.system('rm %s/%s_temp* \n'%(path, gsm))
            gzip_fastq(lay_type,path,gsm)
            gunzip_fastq(lay_type,path,gsm)
        return True
    else:
        return False

# def downloadFastqFromGEO(path,gsm,srr,lay_type):
#     cat_file1 = ''
#     cat_file2 = ''
#     ncbi_link = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR'
#     sys.stdout.write(str(len(srr))+' runs\n')
#     if lay_type == 'SINGLE':
#         for i in range(len(srr)):
#             ftp = ncbi_link + '/' + srr[i][:6] + '/' + srr[i] +'/' + srr[i]+ '.sra'
#             # download the sra files and transmit them into fastq files
#             fsra = '%s/%s_temp%s.sra'%(path, gsm, i+1)
#             os.system('wget %s -O %s \n'%(ftp, fsra))
#             if checkFastqSize(fsra) != True:
#                 return False
#             os.system('echo "+++fastq-dump++++"')
#             os.system('\nfastq-dump %s/%s_temp%s.sra -O %s \n'%(path, gsm,i+1, path))
#             cat_file1 = cat_file1 + '%s/%s_temp%s.fastq '%(path, gsm,i+1) 
#         os.system('cat %s> %s/%s.fastq \n'%(cat_file1,path,gsm))
#         os.system('rm %s/%s_temp* \n'%(path, gsm))
#         gzip_fastq(lay_type,path,gsm)
#     elif lay_type == 'PAIRED':
#         for i in range(len(srr)):
#             ftp = ncbi_link + '/' + srr[i][:6] + '/' + srr[i] +'/' + srr[i]+ '.sra'
#             fsra = '%s/%s_temp%s.sra'%(path, gsm, i+1)
#             os.system('wget %s -O %s \n'%(ftp, fsra))
#             if checkFastqSize(fsra) != True:
#                 return False
#             os.system('echo "+++fastq-dump++++" \n')
#             os.system('\nfastq-dump --split-files %s/%s_temp%s.sra -O %s \n'%(path, gsm,i+1, path))
#             cat_file1 = cat_file1 + '%s/%s_temp%s_1.fastq '%(path, gsm,i+1)
#             cat_file2 = cat_file2 + '%s/%s_temp%s_2.fastq '%(path, gsm,i+1)
#         os.system('cat %s> %s/%s_R1.fastq \n'%(cat_file1, path,gsm))
#         os.system('cat %s> %s/%s_R2.fastq \n'%(cat_file2, path,gsm))
#         os.system('rm %s/%s_temp* \n'%(path, gsm))
#         gzip_fastq(lay_type,path,gsm)
#     return True

def downloadMicroarrayData(gsm_html):
    link_regexp = re.compile(r'href=.*[Cc][Ee][Ll].*>\(ftp\)')
    link_infor = link_regexp.search(gsm_html)
    if link_infor:
        link = link_infor.group().replace('">(ftp)','').replace('href="','')
        sys.stdout.write("Found ftp download link.")
        cmd = "wget %s" % link
        return cmd
    else:
        sys.stderr.write("!!!!! ERROR: Failed to download find CEL files. !!!!!\n")
        sys.exit(6)

def getLayType(srx_html,gsm):
    lay_type = re.search('<div>Layout: <span>.{6}</span>',srx_html)
    if lay_type:
        lay_type = lay_type.group()
        lay_type = lay_type[-13:-7]
        sys.stderr.write(lay_type + "\n")
        if lay_type != "SINGLE" and lay_type != "PAIRED":
            lay_type = None
            sys.stderr.write('!!!!! ERROR: neither PAIRED nor SINGLE end sequencing: %s. !!!!!\n'%gsm)
            sys.exit(3)
        else:
            return lay_type

def getSRR(srx_html):
    srr_regexp = re.compile('>SRR[0-9]*</a></td><td')
    srr = srr_regexp.findall(srx_html)
    srr = [i.lstrip(">").split("</a")[0] for i in srr]
    return srr

def getSrxHtml(gsm_html):
    srx_infor = re.search('https://www.ncbi.nlm.nih.gov/sra\S*"',gsm_html)
    if srx_infor:
        srx = srx_infor.group().rstrip('"').lstrip('https://www.ncbi.nlm.nih.gov/sra?term=')
        # get the SRR id('>SRR1588518</a></td><td') and find the type of layout
        srx_url = 'http://www.ncbi.nlm.nih.gov/sra?term=%s'%srx
        srx_html = urllib.request.urlopen(srx_url).read().decode('utf-8')
        return srx_html
    else:
        return None

def getGsmHtml(gsm):
    gsm_url = 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s'%gsm
    gsm_html = urllib.request.urlopen(gsm_url).read().decode('utf-8')
    return gsm_html


def main():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser()
    parser.add_argument('-i', '--id', help='GSM id', required=True)
    parser.add_argument('-o', '--output', help='the path to save fastq', default=".")
    parser.add_argument('-g', '--gzip', help='the flag of whether compress fastq files', action = "store_true", default=False)
    args = parser.parse_args()

    gsm = args.id
    path = args.output
    global compress
    compress = args.gzip

    if compress == True:
        sys.stdout.write("Gzip files.\n")
    else:
        sys.stdout.write("Would not compress files.\n")

    try:
        os.makedirs(path)
        sys.stdout.write("Folder %s was created!\n" % path)
    except FileExistsError:
        sys.stdout.write("Folder %s existed!\n" % path)

    os.system('echo %s' % gsm)
    gsm_html = getGsmHtml(gsm)
    srx_html = getSrxHtml(gsm_html)
    if srx_html:
        SRR = getSRR(srx_html)
        if SRR:
            LayType = getLayType(srx_html,gsm)
            # GEO_status = downloadFastqFromGEO(path,gsm,SRR,LayType)
            # if GEO_status == False:
            #     sys.stderr.write("+++++ WARNING: Could NOT download %s sra file from GEO FTP server. Trying Prefetch (sra-tools). +++++\n\n" % gsm)
            PREFETCH_status = downloadFastqByPrefetch(path,gsm,SRR,LayType)
            if PREFETCH_status == False:
                sys.stderr.write("+++++ WARNING: Could NOT download %s sra file by PREFETCH. Trying EBI. +++++\n\n" % gsm)
                # EBI_status = downloadFastqFromEBI(path,gsm,SRR,LayType)
                # if EBI_status == False:
                #     sys.stderr.write("!!!!! ERROR: Could NOT download %s sra file from EBI. END QUERY. !!!!!\n\n" % gsm)
                sys.exit(5)
        else:
            sys.stderr.write("!!!!! ERROR: Do not find SRR information for %s. Experiment data may not be public. !!!!!\n\n" % gsm)
            sys.exit(4)
    else: # may not seq data
        sys.stderr.write("+++++ WARNING: Do not find SRX information for %s. Microarray Data? +++++\n\n" % gsm)
        cmd = downloadMicroarrayData(gsm_html)
        cmd += " -P %s" % path
        os.system(cmd)


if __name__ == "__main__":
    main()
    






