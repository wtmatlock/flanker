from flanker.flanker import *
from flanker.cluster import *



def flank_salami_linear(file, include_gene,step, stop,gene,flank):
    args = get_arguments()
    unfiltered_abricate_file = str(file + '_resfinder') # name of abricate output for fasta
    data = pd.read_csv(unfiltered_abricate_file, sep='\t', header = 0)

    guids=data['SEQUENCE'].unique()

    for guid in guids:
        abricate_file=filter_abricate(data,guid)
        pos = flank_positions(abricate_file, gene)
        if pos == True:
            log.error(f"Error: Gene {gene} not found in {guid}")


        else:
             gene_sense=abricate_file.loc[abricate_file['GENE']==gene].filter(items=['STRAND'])



             gene_sense=str(gene_sense['STRAND'].iloc[0])






             d_lin = {(True, 'upstream'): lambda record, positions, w, l : record.seq[max(0,positions[0]):min(l,positions[1])],
             (True, 'downstream'): lambda record, positions, w, l : record.seq[positions[0]:min(l, positions[1])],
             (False, 'upstream'): lambda record, positions, w, l : record.seq[max(0, positions[0]):positions[0]],
             (False, 'downstream'): lambda record, positions, w, l : record.seq[positions[1]:min(l, positions[1]+w)]}

             w = int(step)
             x = args.flank
             start_left=pos[0]
             start_right=pos[1]

             print(start_left,start_right,step)

             if include_gene == True:
                 start_left=pos[1]
                 start_right=pos[0]




             for record in SeqIO.parse(file, "fasta"):


                 if record.description == guid:

                     if gene_sense == '-':

                         #record.seq = record.seq.reverse_complement()
                         if flank == 'upstream':
                             x = 'downstream'
                         else:
                             x = 'upstream'
                     name=str(record.description)

                     log.info(f"{gene} found in {record.description}")

                     l = len(record.seq)





                     for i in range(0,stop,step):
                         print(start_left)
                         print(start_right)

                         for record in SeqIO.parse(file, "fasta"):
                             name=str(record.description)
                             print(pos[2] + ' found')

                             s = int(start)
                             l = len(record.seq)

                             if flank == 'upstream':


                                 record.seq = record.seq[max(0, start_left-step):start_left]

                                 record.description = f"{record.description} | {pos[2]} | {step}bp window"

                                 with open(f"{name}_{pos[2]}_{i}_salami_left_flank.fasta", "w") as f:
                                     SeqIO.write(record, f, "fasta")
                                     print(f"{f.name} sucessfully created!")
                                     f.close()
                                 start_left=start_left-step


                             elif flank == 'downstream':


                                record.seq = record.seq[start_right:min(len(record.seq), start_right + step)]

                                record.description = f"{record.description} | {pos[2]} | {step}bp window"

                                with open(f"{name}_{pos[2]}_{i}_salami_right_flank.fasta", "w") as f:
                                    SeqIO.write(record, f, "fasta")
                                    print(f"{f.name} sucessfully created!")
                                    f.close()
                                start_right=start_right+step





def salami_main(genes,fasta,include_gene,wstep,wstop,out,flank,threads,threshold,cluster):
    with open(genes) as gene_list:
        for gene in gene_list:
            print(gene)
            run_abricate(fasta)

            print("Working on gene {}".format(gene))
            flank_salami_linear(fasta,include_gene,wstep,wstop,gene.strip(),flank)
            print(cluster)
        if cluster ==True:
            define_clusters("salami","mode",threads,threshold,out)
            log.info("Cleaning up")
            flank_scrub()
