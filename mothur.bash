sh
PROJECTNAME=$1




mothur "#make.contigs(file=$PROJECTNAME.file, processors=32);
count.groups(count=current);
summary.seqs(fasta=current, count=current);
screen.seqs(fasta=current, count=current, summary=current, maxambig=0, maxlength=275);
summary.seqs(fasta=current, count=current);
unique.seqs(fasta=current, count=current);
summary.seqs(fasta=current, count=current);
align.seqs(fasta=current, reference=silva.nr_v138_1.v4.align);
summary.seqs(fasta=current, count=current);
screen.seqs(fasta=current, count=current, summary=current, start=2, end=8580, maxhomop=8);
filter.seqs(fasta=current, vertical=T);
summary.seqs(fasta=current, count=current);
pre.cluster(fasta=current, diffs=2, count=current);
summary.seqs(fasta=current, count=current);
chimera.vsearch(fasta=current, count=current, dereplicate=t);
count.groups(count=current);
summary.seqs(fasta=current, count=current);
classify.seqs(fasta=current, count=current, reference=silva.nr_v138_1.v4.align, taxonomy=silva.nr_v138_1.tax, cutoff=80);
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Eukaryota);
summary.tax(taxonomy=current, count=current);
dist.seqs(fasta=current, countends=F, cutoff= 0.03, processors=16);
cluster(column=current, count=current, method=opti);
summary.seqs(processors=32);
make.shared(list=current, count=current);
classify.otu(list=current, count=current, taxonomy=current);
get.oturep(fasta=current, count=current, list=current, method=abundance);
count.groups(shared=current);
rarefaction.single(shared=current);
summary.single(shared=current, calc=nseqs-sobs-coverage-shannon-shannoneven-invsimpson, subsample=10000);
dist.shared(shared=current, calc=braycurtis-jclass-thetayc, subsample=10000);
sub.sample(count=$PROJECTNAME.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.count_table, shared=current, list=$PROJECTNAME.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.opti_mcc.list, size=10000, persample=true, label=0.03);

sub.sample(taxonomy=$PROJECTNAME.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.nr_v138_1.wang.pick.taxonomy, count=$PROJECTNAME.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.count_table, list=$PROJECTNAME.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.opti_mcc.list, size=10000, persample=true, label=0.03);
summary.tax(taxonomy=current, count=current); system(mkdir send);
system(cp *shared send); system(cp *cons.tax* send); system(cp *pick.tax.summary send); system(cp *pick.subsample.tax.summary send); system(cp *.rep.fasta send); system(cp *lt.ave.dist send); system(cp *groups.ave-std.summary send); system(cp mothur.bash send); system(cp mothur.*.logfile send);"

