from flanker.flanker import *
from flanker.cluster import *

def flank_salami_linear(start,stop,step,file,gene):
        args = get_arguments()
        abricate_file = str(file + '_resfinder') # name of abricate output for fasta

        pos = flank_positions(abricate_file, gene)

        if pos != True:
            start=pos[0]
            start_right=pos[1]
            for i in range(0,stop,step):
                print(i)

                for record in SeqIO.parse(file, "fasta"):
                    name=str(record.description)
                    print(pos[2] + ' found')

                    s = int(start)
                    l = len(record.seq)

                    if args.flank == 'left':


                        record.seq = record.seq[max(0, start-step):start]

                        record.description = f"{record.description} | {pos[2]} | {step}bp window"

                        with open(f"{name}_{pos[2]}_{i}_salami_left_flank.fasta", "w") as f:
                            SeqIO.write(record, f, "fasta")
                            print(f"{f.name} sucessfully created!")
                            f.close()



                            start=start-step

                    elif args.flank == 'right':


                        record.seq = record.seq[start_right:min(len(record.seq), start_right + step)]

                        record.description = f"{record.description} | {pos[2]} | {step}bp window"

                        with open(f"{name}_{pos[2]}_{i}_salami_right_flank.fasta", "w") as f:
                            SeqIO.write(record, f, "fasta")
                            print(f"{f.name} sucessfully created!")
                            f.close()



                            start_right=start_right+step

def salami_main(genes,fasta,window,wstep,wstop,indir,out,threads,threshold,cluster):
    with open(genes) as gene_list:
        for gene in gene_list:

            run_abricate(fasta)

            print("Working on gene {}".format(gene))
            flank_salami_linear(window,wstop,wstep,fasta,gene.strip())
            print(cluster)
        if cluster ==True:
            define_clusters(gene,"Salami",indir,threads,threshold,out)

            filelist=glob.glob(str(indir + str("*flank.fasta")))
            #print("Cleaning")
            #for filename in filelist:
                #os.remove(filename)
                #print("Here we go again")
