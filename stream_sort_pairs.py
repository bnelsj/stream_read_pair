from optparse import OptionParser
import pysam
import pdb
from sys import stdout
import subprocess
import numpy as np
import random
import pdb 

def fetch_all(b):
    for l in b.fetch(until_eof=True):
        yield l

def sam_str(r):
    
    return "%s\n"%"\t".join([r.qname,
                     str(r.flag),
                     str(r.rname+1),
                     str(r.pos+1),
                     str(r.mapq),
                     r.cigarstring,
                     r.rnext==-1 and "*" or str(r.rnext),
                     r.rnext==-1 and "0" or str(r.pnext),
                     str(r.tlen),
                     r.seq,
                     str(r.qual)]+[":".join([str(t) for t in tag]) for tag in r.tags]+["QS:%d"%(r.qstart),"QE:%d"%(r.qend)])
                     
class pairing_window(object):

    def __init__(self, wnd_size=100000):
        self.curr_contig = None
        self.wnd_size = wnd_size
        self.wnd_start = None
        self.wnd_end = None
        self.reads_by_pos = {}
        self.reads_by_name = {}
        self.n_pairs_output = 0
    
    def clean_up_all(self):
        del self.reads_by_pos
        del self.reads_by_name
        self.reads_by_pos = {}
        self.reads_by_name = {}
        
    def update_wnd(self, e_pos):
        self.wnd_end = e_pos
        while self.wnd_end - self.wnd_start >self.wnd_size:
            if self.wnd_start in self.reads_by_pos:
                for read in self.reads_by_pos[self.wnd_start]:
                    del self.reads_by_name[read.qname]
                del self.reads_by_pos[self.wnd_start]
            self.wnd_start+=1
        
    def add_read(self, read, binary):
        """
        reset
        """
        if read.rname != self.curr_contig:
            self.clean_up_all()
            self.curr_contig = read.rname
            self.wnd_start = read.pos
            self.wnd_end = read.pos+self.wnd_size
        
        if read.qname in self.reads_by_name:
            """
            ouput read
            """
            if binary:
                outstream.write(self.reads_by_name[read.qname])
                outstream.write(read)
            else:
                outstream.write(sam_str(self.reads_by_name[read.qname]) + '\n')
                outstream.write(sam_str(read) + '\n')
            self.n_pairs_output +=1
        else:
            """
            add read 
            """
            self.reads_by_name[read.qname] = read
            if not read.pos in self.reads_by_pos:
                self.reads_by_pos[read.pos] = []
            self.reads_by_pos[read.pos].append(read)
        
        if read.pos > self.wnd_end:
            self.update_wnd(read.pos)
       
def is_good_read(read, bamfile, contigs_to_consider):
    if not read.is_secondary and not read.is_qcfail and not read.is_duplicate and not read.is_unmapped and bamfile.getrname(read.tid) in contigs_to_consider:
        return True
    else:
        return False
 
if __name__=="__main__":

    opts = OptionParser()
    opts.add_option('','--input_bam',dest='fn_bam', default='')
    opts.add_option('','--include_chrs', default='', help='list of chrs to consider separated by :')
    opts.add_option('','--window',dest='window', default=100000, type = int)
    opts.add_option('','--n_samples',dest='n_samples', default=100000, type = int)
    opts.add_option('','--subsample_reads',dest='subsample_reads', default=False, action="store_true")
    opts.add_option('','--binary', action='store_true', default=False, help='Write to stream in bam format')

    (o, args) = opts.parse_args()

    if o.fn_bam == '':
        print('Must specify input bam')
        sys.exit(1)
    if o.include_chrs == '':
        print('Must specify chrs to consider')
        sys.exit(1)

    b = pysam.Samfile(o.fn_bam, 'rb')

    if o.binary:
        outstream = pysam.Samfile('/dev/stdout', 'wb', template=b)
    else:
        outstream = open('/dev/stdout', 'w')

    contigs_to_consider = o.include_chrs.split(':')

    if not any(map(lambda x: x.startswith('chr'), b.references)):
        contigs_to_consider = map(lambda x: x.replace('chr', ''), contigs_to_consider)

    if o.subsample_reads:
        contigs_to_len = {contig:int(b.lengths[i]) for i, contig in enumerate(b.references) if contig in contigs_to_consider}
        contigs_to_start = {}
        
        for contig in contigs_to_consider:
            read = None
            for r in b.fetch(contig):
                read = r
                break
            if read == None:
                del contigs_to_len[contig]
            else:
                contigs_to_start[contig] = read.pos

        t = np.sum(np.array(contigs_to_len.values()))
      
        for contig, l in contigs_to_len.iteritems():
            n_reads = int(o.n_samples * (l / float(t)))
            ps = sorted([random.randrange(contigs_to_start[contig], l) for i in xrange(n_reads)])
            
            for p in ps:
                pairing_obj = pairing_window(wnd_size = o.window) 
                curr_c = pairing_obj.n_pairs_output

                for read in b.fetch(contig, p, l):
                    if is_good_read(read, b, contigs_to_consider):
                        pairing_obj.add_read(read, o.binary)
                        if pairing_obj.n_pairs_output > curr_c + 100:
                            break
                del pairing_obj
    else:
        pairing_obj = pairing_window(wnd_size = o.window) 
        for read in fetch_all(b):
            if is_good_read(read, b, contigs_to_consider):
                pairing_obj.add_read(read, o.binary)

    b.close()
    outstream.close()
