import os, sys
import urllib.request, urllib.error, urllib.parse
import re
import argparse

# transfer the sra files to fastqfiles and cat the multiple fastq files
def catfastq(main_path, gsm,srx_infor,srr,lay_type):
    cmd = ''
    cat_file1 = ''
    cat_file2 = ''
    srx_infor = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR'
    print((lay_type+'\n'+str(len(srr))+' runs'))
    if lay_type == 'SINGLE':
        for i in range(len(srr)):
            srr[i] = srr[i][1:-12]
            #ftp = srx_infor + '//' + srr[i] + '/' + srr[i] + '.sra'
            ftp = srx_infor + '/' + srr[i][:6] + '/' + srr[i] +'/' + srr[i]+ '.sra'
            # download the sra files and transmit them into fastq files
            fsra = '%s/%s_temp%s.sra'%(main_path, gsm, i+1)
            cmd = cmd + 'wget %s -O %s \n'%(ftp, fsra) 
            cmd = cmd + 'echo "+++fastq-dump++++" \n\nfastq-dump %s/%s_temp%s.sra -O %s \n'%(main_path, gsm,i+1, main_path) 
            cat_file1 = cat_file1 + '%s/%s_temp%s.fastq '%(main_path, gsm,i+1) 
        cmd = cmd + 'cat %s> %s/%s.fastq \n'%(cat_file1,main_path,gsm) 
        cmd = cmd + 'rm %s/%s_temp* \n'%(main_path, gsm)
        cmd = cmd + 'gzip %s/%s.fastq \n' % (main_path, gsm)
    elif lay_type == 'PAIRED':
        for i in range(len(srr)):
            srr[i] = srr[i][1:-12]
            ftp = srx_infor + '/' + srr[i][:6] + '/' + srr[i] +'/' + srr[i]+ '.sra'
            fsra = '%s/%s_temp%s.sra'%(main_path, gsm, i+1)
            cmd = cmd + 'wget %s -O %s \n'%(ftp, fsra) 
            # cmd = cmd + 'python '+os.path.join(main_path+'/modules/checkSRA.py')+' %s\nif [ $? == 1 ];then\necho %s >> %s/SRAchecking_fail.xls\nexit\nfi\n'%(fsra, gsm+','+ftp+','+fsra, recordFolder)
            cmd = cmd + 'echo "+++fastq-dump++++" \n\nfastq-dump --split-files %s/%s_temp%s.sra -O %s \n'%(main_path, gsm,i+1, main_path) 
            cat_file1 = cat_file1 + '%s/%s_temp%s_1.fastq '%(main_path, gsm,i+1)
            cat_file2 = cat_file2 + '%s/%s_temp%s_2.fastq '%(main_path, gsm,i+1)
        cmd = cmd + 'cat %s> %s/%s_R1.fastq \n'%(cat_file1, main_path,gsm) 
        cmd = cmd + 'cat %s> %s/%s_R2.fastq \n'%(cat_file2, main_path,gsm) 
        cmd = cmd + 'rm %s/%s_temp* \n'%(main_path, gsm)
        cmd = cmd + 'gzip %s/%s_R1.fastq \n'%(main_path, gsm)
        cmd = cmd + 'gzip %s/%s_R2.fastq \n'%(main_path, gsm)
    return(cmd)



def LinkPlusDownload(gsm, path):
    """get link sequece type and return the path of sbatch file
    """
    os.system('echo %s'%gsm)
    gsm_url = 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s'%gsm
    gsm_handler = urllib.request.urlopen(gsm_url)
    gsm_html = gsm_handler.read()
    gsm_html = gsm_html.decode('utf-8')
    # get the ftp location of SRX file and the SRX id
    srx_regexp = re.compile('https://www.ncbi.nlm.nih.gov/sra\S*"')
    srx_infor = srx_regexp.search(gsm_html)
    if srx_infor:
        srx = srx_infor.group().rstrip('"').lstrip('https://www.ncbi.nlm.nih.gov/sra?term=')
    else:
        print(('%s: no srx number'%gsm))
    # get the SRR id('>SRR1588518</a></td><td') and find the type of layout
    srx_url = 'http://www.ncbi.nlm.nih.gov/sra?term=%s'%srx
    srx_handler = urllib.request.urlopen(srx_url)
    srx_html = srx_handler.read()
    srx_html = srx_html.decode('utf-8')
    # find the layout type (<div>Layout: <span>SINGLE</span>)
    lay_infor = re.compile('<div>Layout: <span>.{6}</span>')
    lay_type = lay_infor.search(srx_html)
    lay_type = lay_type.group()
    lay_type = lay_type[-13:-7]
    # get the srr id and download the sra files
    srr_regexp = re.compile('>SRR[0-9]*</a></td><td')
    srr = srr_regexp.findall(srx_html)
    if lay_type == 'SINGLE':
        fastq_file = '%s/%s.fastq'%(path,gsm)
         # single.append(gsm)
         # print >>SingleEnd_gsm, gsm
    elif lay_type == 'PAIRED':
        fastq_file = '%s/%s_R1.fastq,%s/%s_R2.fastq'%(path,gsm,path,gsm)
    else:
        print(('+++neither PAIRED nor SINGLE end sequencing: %s+++'%gsm))
        return(None, None)
    fastCmd = catfastq(path, gsm,srx_infor,srr,lay_type)
    os.system(fastCmd)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("++download fastq data by GSM id:", file=sys.stderr)
        print("gsmid\tpath", file=sys.stderr)
        print("1\tgsmid", file=sys.stderr)
        print("2\toptional: the path to save fastq", file=sys.stderr)
        sys.exit(1)
    elif len(sys.argv) == 2:
    	gsm = sys.argv[1]
    	path = os.getcwd()
    else:
    	gsm = sys.argv[1]
    	path = sys.argv[2]
    LinkPlusDownload(gsm, path)




