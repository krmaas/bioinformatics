sh
PROJECTNAME=$1





mothur "#make.contigs(file=$PROJECTNAME.file, processors=32);
summary.seqs(fasta=current, count=current);
count.groups(count=current);
screen.seqs(fasta=current, count=current, summary=current, maxambig=2, maxlength=400);
summary.seqs(fasta=current, count=current);
unique.seqs(fasta=current, count=current);
summary.seqs(fasta=current, count=current);
pre.cluster(fasta=current, diffs=2, count=current);
summary.seqs(fasta=current, count=current);
chimera.vsearch(fasta=current, count=current, dereplicate=t);
count.groups(count=current);
summary.seqs(fasta=current, count=current);
classify.seqs(fasta=current, count=current, reference=UNITEv6_sh_dynamic.fasta, taxonomy=UNITEv6_sh_dynamic.tax, cutoff=60);
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=unknown-Protista);
summary.tax(taxonomy=current, count=current);
cluster(fasta=current, count=current, method=agc, cutoff=0.05);
summary.seqs(processors=32);
make.shared(list=current, count=current);
classify.otu(list=current, count=current, taxonomy=current);
get.oturep(fasta=current, count=current, list=current, method=abundance);
count.groups(shared=current);
rarefaction.single(shared=current);
summary.single(shared=current, calc=nseqs-sobs-coverage-shannon-shannoneven-invsimpson, subsample=10000);
dist.shared(shared=current, calc=braycurtis-jclass-thetayc, subsample=10000);
sub.sample(count=$PROJECTNAME.trim.contigs.good.unique.precluster.denovo.vsearch.pick.count_table, shared=current, list=$PROJECTNAME.trim.contigs.good.unique.precluster.denovo.vsearch.pick.agc.list, size=10000, persample=true, label=0.05);
sub.sample(taxonomy=$PROJECTNAME.trim.contigs.good.unique.precluster.denovo.vsearch.UNITEv6_sh_dynamic.wang.pick.taxonomy, count=$PROJECTNAME.trim.contigs.good.unique.precluster.denovo.vsearch.pick.count_table, list=$PROJECTNAME.trim.contigs.good.unique.precluster.denovo.vsearch.pick.agc.list, size=10000, persample=true, label=0.05);
summary.tax(taxonomy=current, count=current); system(mkdir send);
system(cp *shared send); system(cp *cons.tax* send); system(cp *pick.tax.summary send); system(cp *pick.subsample.tax.summary send); system(cp *.rep.fasta send); system(cp *lt.ave.dist send); system(cp *groups.ave-std.summary send); system(cp mothur.bash send); system(cp mothur.*.logfile send);"
