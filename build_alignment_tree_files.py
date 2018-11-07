import os
from Bio import SeqIO
from ete3 import Tree


def build_files(input_dir, gene):
    print('.' * 10 + gene + '.' * 10)
    gene_dir = os.path.join(input_dir, gene)
    cds_file = os.path.join(gene_dir, '{0}.cds'.format(gene))
    out_group = ''
    for seq_record in SeqIO.parse(cds_file, 'fasta'):
        if 'Atri' in seq_record.id:
            out_group = seq_record.id
            break
    faa_file = os.path.join(gene_dir, '{0}.fasta'.format(gene))
    faa_aln_file = os.path.join(gene_dir, '{0}_aln.fasta'.format(gene))
    mafft_cmd = 'mafft --auto --quiet {0} > {1}'.format(faa_file, faa_aln_file)
    os.system(mafft_cmd)
    cds_aln_file = os.path.join(gene_dir, '{0}_aln.paml'.format(gene))
    cds_aln_cmd = 'perl {0} {1} {2} -output paml > {3}'.format(pal2nal,
                                                               faa_aln_file,
                                                               cds_file,
                                                               cds_aln_file)
    os.system(cds_aln_cmd)
    os.chdir(gene_dir)
    raxml_cmd = 'raxmlHPC-PTHREADS-SSE3 -m PROTGAMMAAUTO -p 12345 -x 12345 ' \
                '-# 500 -T 18 -s {0} -n nwk -f a ­­silent'.format(faa_aln_file)
    os.system(raxml_cmd)
    tmp_tree = os.path.join(gene_dir, 'RAxML_bestTree.nwk')
    rooted_tree = os.path.join(gene_dir, 'rooted_{0}.nwk'.format(gene))
    t = Tree(tmp_tree, format=1)
    t.set_outgroup(out_group)
    t.write(outfile=rooted_tree, format=1)
    os.chdir(my_path)


if __name__ == '__main__':
    my_path = os.getcwd()
    pal2nal = os.path.join(my_path, 'pal2nal.pl')
    all_genes_dir = os.path.join(my_path, 'all_plant_genes')
    for root, dirs, files in os.walk(all_genes_dir):
        for each_dir in dirs:
            build_files(all_genes_dir, each_dir)
