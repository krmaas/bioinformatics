#!bash
PROJECTNAME=$1

mkdir -p $PROJECTNAME/fastq
for f in basespace/Projects/$PROJECTNAME/Samples/*/Files/*.gz; 
do s=${f##basespace/Projects/$PROJECTNAME/Samples/}; s=${s%%/*}; 
cp $f $PROJECTNAME"/fastq/"$s"."${f##*Files/}; 
done

cd $PROJECTNAME

mv fastq/fileList.paired.file $PROJECTNAME.file

module load mothur/1.38.1
mothur "#make.file(inputdir=fastq, type=gz);"

module load mothur/1.39.4

mothur "#make.contigs(file=$PROJECTNAME.file); make.contigs(processors=32, file=current); summary.seqs(fasta=current); screen.seqs(fasta=current, group=current, summary=current, maxambig=0, maxlength=275); summary.seqs(fasta=current); unique.seqs(fasta=current); summary.seqs(fasta=current, name=current); count.seqs(name=current, group=current); align.seqs(fasta=current, reference=silva.nr_v119.v4.align); summary.seqs(fasta=current, count=current); screen.seqs(fasta=current, count=current, summary=current, start=1968, end=11550, maxhomop=8); filter.seqs(fasta=current, vertical=T); summary.seqs(fasta=current, count=current);pre.cluster(fasta=current, diffs=2, count=current); summary.seqs(fasta=current, count=current); chimera.uchime(fasta=current, count=current, dereplicate=t); remove.seqs(fasta=current, accnos=current, count=current); summary.seqs(fasta=current, count=current); classify.seqs(fasta=current, count=current, reference=silva.nr_v119.v4.align, taxonomy=silva.nr_v119.tax, cutoff=80); remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota); summary.tax(taxonomy=current, count=current); dist.seqs(fasta=current, countends=F, cutoff= 0.03, processors=16); cluster(column=current, count=current, method=opti); summary.seqs(processors=32); make.shared(list=current, count=current); classify.otu(list=current, count=current, taxonomy=current); get.oturep(fasta=current, count=current, list=current, method=abundance); count.groups(shared=current); summary.single(shared=current, calc=nseqs-sobs-coverage-shannon-shannoneven-invsimpson, subsample=10000); dist.shared(shared=current, calc=braycurtis-jest-thetayc, subsample=10000); sub.sample(taxonomy=$PROJECTNAME.trim.contigs.good.unique.good.filter.precluster.pick.nr_v119.wang.pick.taxonomy, count=$PROJECTNAME.trim.contigs.good.unique.good.filter.precluster.denovo.uchime.pick.pick.pick.count_table, list=$PROJECTNAME.trim.contigs.good.unique.good.filter.precluster.pick.pick.opti_mcc.list, size=10000, persample=true, label=0.03); summary.tax(taxonomy=current, count=current); system(mkdir send); system(cp *shared send); system(cp *cons.tax* send); system(cp *pick.tax.summary send); system(cp *pick.subsample.tax.summary send); system(cp *.rep.fasta send); system(cp *lt.ave.dist send); system(cp *groups.ave-std.summary send); system(cp mothur.batch send); system(cp mothur.*.logfile send);"

