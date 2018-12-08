from Bio import SeqIO
import itertools 
import pandas as pd
import glob
import networkx as nx


SDA = "/net/eichler/vol21/projects/bac_assembly/nobackups/genomeWide/CHM13_V2/LocalAssemblies/SDA.contigs.fasta"
PSV = "/net/eichler/vol21/projects/bac_assembly/nobackups/genomeWide/CHM13_V2/LocalAssemblies/psv.tbl"
#READS = "/net/eichler/vol2/home/mvollger/ontultralong/ont.ul.fofn"
#reads = "/net/eichler/vol21/projects/bac_assembly/nobackups/genomeWide/CHM13/LocalAssemblies/ont.200.kbp.fastq.fofn"
READS = "ont.ul.fofn"

minscore = 10000
minlenread = 150000


rule all:
	input:
		final="final",


rule MapSDAtoRef:
	input:
		ref="/net/eichler/vol2/home/mvollger/assemblies/hg38/ucsc.hg38.no_alts.fasta",
		reads=SDA,
	output:
		bam="PlaceOnt/SDAtoRef.bam"
	threads:6
	shell:"""
/net/eichler/vol2/home/mvollger/projects/minimap2/minimap2 -t {threads} -2 --eqx -L -r 50000 -ax asm10 \
		-R '@RG\\tID:minimap\\tSM:FAKE_SAMPLE\\tPL:ont' \
		{input.ref} {input.reads} | \
		samtools view -@ {threads} -bS -F 2308 - | \
		samtools sort -@ {threads} -m 8G -T tmp -o {output.bam}
		samtools index {output.bam}
"""
onts = open(READS).readlines()

IDS = list(range(len(onts)))
ONTS = {}
for idx, ont in enumerate(onts):
	ONTS[idx] = ont.strip()
#print(ONTS, IDS)

def ontfiles(ID):
	ID = int(str(ID))
	return(ONTS[ID])


rule GetLargeONT:
	input:
		fastq = ontfiles, 
	output:
		fasta = "LargeONT/{ID}.fasta",
	run:
		recs = list(SeqIO.parse(input["fastq"], "fastq"))
		IDs = []
		keep = []
		t = len(recs)
		for idx, rec in enumerate(recs):
			if( (rec.id not in IDs) and (len(rec.seq) > minlenread) ):
				IDs.append(rec.id)
				keep.append(rec)
			if(idx%10000 == 0):
				print(idx, t)
		SeqIO.write(keep, output["fasta"], "fasta")

rule MapOntToRef:
	input:
		ref="/net/eichler/vol2/home/mvollger/assemblies/hg38/ucsc.hg38.no_alts.fasta",
		reads=expand("LargeONT/{ID}.fasta", ID=IDS),
	output:
		bam="PlaceOnt/ontToRef.bam"
	threads:20
	shell:"""
/net/eichler/vol2/home/mvollger/projects/minimap2/minimap2 -t {threads} -2 --eqx -L -r 50000 -ax map-ont \
		-R '@RG\\tID:minimap\\tSM:FAKE_SAMPLE\\tPL:ont' \
		{input.ref} {input.reads} | \
		samtools view -@ {threads} -bS -F 2308 - | \
		samtools sort -@ {threads} -m 8G -T tmp -o {output.bam}
		samtools index {output.bam}
"""




rule MapONTtoSDA:
	input:
		sda = SDA,
		bam="PlaceOnt/ontToRef.bam",
		reads=expand("LargeONT/{ID}.fasta", ID=IDS),
	output:
		sdaToOnt = "PlaceOnt/OntToSDA.nof.bam",
	threads:12
	shell:"""
minimap2 -ax map-ont -L -2 -t {threads} -r 5000 -s {minscore} -N 30 -p .10 --eqx \
	{input.sda} {input.reads} > tmp.sam

samtools view -b -F 2052 tmp.sam | \
	samtools sort -m 4G -T tmp -o {output.sdaToOnt}
	samtools index {output.sdaToOnt}
"""

rule MapSDAtoSDA:
	input:
		sda = SDA,
	output:
		sdabam = "PlaceOnt/SDAtoSDA.bam"
	threads:12
	shell:"""
minimap2 -t {threads} -r 50000 -ax asm20 -s 8000 -N 10 -p .05 --eqx \
		{input.sda} {input.sda} | \
		samtools view -b -F 4 - | \
		samtools sort -m 4G -T tmp -o {output.sdabam}
		samtools index {output.sdabam}
"""

rule filterAlns:
	input:
		ontToSda= "PlaceOnt/OntToSDA.nof.bam",
	output:
		filt = "PlaceOnt/OntToSDA.filter.bam",
		lhr = "PlaceOnt/LHR.tsv",
	shell:"""
~/projects/SDA/PSVFilterAlignments.py --minPSVs 5 --minFRAC 0.75 --minLHR 0 --bam {input.ontToSda} --psvs psv.tbl --values {output.lhr} {output.filt}
"""

IDS = ["nof", "filter"]

rule makeTBL:
	input:
		ontToSda= "PlaceOnt/OntToSDA.{ID}.bam",
	output:
		tsv = "PlaceOnt/OntToSDA.{ID}.bam.tsv",
	shell:"""
~/projects/SDA/scripts/samIdentity.py --header {input.ontToSda} > {output.tsv}
"""

rule MakeBed:
	input:
		ont = "PlaceOnt/ontToRef.bam",
		sda = "PlaceOnt/SDAtoRef.bam",
	output:
		ont = "PlaceOnt/ontToRef.bed",
		sda = "PlaceOnt/SDAtoRef.bed",
	shell:"""
bedtools bamtobed -i {input.ont} > {output.ont}
bedtools bamtobed -i {input.sda} > {output.sda}
"""



rule refintersect:
	input:
		sda = "PlaceOnt/SDAtoRef.bed",
		ont = "PlaceOnt/ontToRef.bed",
	output:
		bed = "PlaceOnt/intersect.bed"	
	shell:"""
# intersect and require 25 kbp of overlap 
bedtools intersect -a {input.ont} -b {input.sda} -wo | \
	awk '$13 > 5000 {{print $0}}' > {output.bed}
"""	



def makeEdges(df2, ont, collapse, realset=None):
		groups = df2.groupby(ont)
		Edges = []
		for name, group in groups:
			if(len(group) != 1):
				# get all pairwise edges to get combinations
				for pair in itertools.combinations(group[collapse], 2):
					pair1 = (pair[0], pair[1], name)
					pair2 = (pair[1], pair[0], name)
					Edges.append( pair1 )
					Edges.append( pair2 )
					#print(row.collapse, row.reference_name)
				#print("")
		Edges = pd.DataFrame(Edges, columns=["node1", "node2", "by"])
		Edges.drop_duplicates(inplace=True)
		return(Edges)


def makeSets(f):
	sets = set()
	lines = open(f).readlines()
	for line in lines:
		if("node1" in line):
			continue 
		token = line.strip().split()
		name1 = token[0] +"_"+ token[1]
		name2 = token[1] +"_"+ token[0]
		sets.add(name1)
		sets.add(name2)
	return(set(sets))

rule alnAcc:
	input:
		bed = "PlaceOnt/intersect.bed",
		bam = "PlaceOnt/OntToSDA.{ID}.bam",
	output:
		acc = "PlaceOnt/AlnCorrect.{ID}.txt",
		sda = temp("PlaceOnt/tmp.sda.{ID}.pairs"),
		ref = temp("PlaceOnt/tmp.ref.{ID}.pairs"),
	run:
		shell("samtools view {input.bam} | cut -f1,3 > {output.sda}")
		shell("cat {input.bed} | cut -f4,10 > {output.ref}")
		sda = makeSets(output["sda"])
		ref = makeSets(output["ref"])
		cor = len(sda.intersection(ref))/2
		total = len(sda)/2
		acc = cor*100.0/total

		out = "Acc:{}\tCorrect:{}\tTotal:{}\n".format(acc, cor, total)
		print(out)
		open(output["acc"], "w+").write(out)

rule refGraph:
	input:
		bed = "PlaceOnt/intersect.bed",
	output:
		edges = "PlaceOnt/ref.edges",
	run:
		df = pd.read_table(input["bed"], header=None)
		Edges = makeEdges(df, 3, 9)	
		Edges.to_csv(output["edges"], sep="\t", index=False)


rule SDAgraph:
	input:
		tsv = "PlaceOnt/OntToSDA.{ID}.bam.tsv",
		edges = "PlaceOnt/ref.edges",
	output:
		sdaedges="PlaceOnt/SDA.{ID}.edges",
	run:
		df = pd.read_table(input["tsv"], sep="\t")
		# only consider links between different ref and qname assicoations, can do by droping duplicaites
		df["collapse"] = df.reference_name.str.extract(r'.*collapse\.(.*)_id\.*', expand=False)
		# sort by collapse and tehn score which says hwo good the alignemtn is. 
		df["score"] = df.matches * 2 - 4*(df.missmatches + df.insertion_events + df.deletion_events)
		df = df.sort_values(by=['reference_name', 'collapse', 'score'])
		# drop pairings for lower scoring alignments from the same collase, keep last only keeps highest score 
		#print(df)

		#df2 = df.drop_duplicates(subset = ["collapse", "reference_name"], keep='last')
		#df2 = df.drop_duplicates(subset = ["reference_name", "collapse"], keep='first')
		
		# read in real edges 
		Edges = makeEdges(df, "query_name", "reference_name")
		Edges.to_csv(output["sdaedges"], sep="\t", index=False)



def makeGraph(df, realset=None):
	nodes = set( list(df.node1) + list(df.node2) )
	g=nx.Graph()
	for nidx, nlabel in enumerate(nodes):
		# add a ndoe with a label
		g.add_node(nlabel, label=nlabel, weight=5)

	for idx, row in df.iterrows():
		if(realset is not None):
			name = row.node1 +"_"+ row.node2
			color="#ff0000"  # red
			weight = 3
			if(name in realset):
				color="#000000" # balck
			g.add_edge(row.node1, row.node2, label=row.by, fill=color, color=color, weight=weight)
		else:
			g.add_edge(row.node1, row.node2, label=row.by)


	# color by connected component
	comps = sorted(nx.connected_components(g), key = len, reverse=True)
	#nx.set_node_attributes(g, "color", "0"*len(g.nodes) )
	outstr = ""
	for idx, cut in enumerate(comps):
		for node in cut:
			g.nodes[node]["color"] = idx
			outstr += g.nodes[node]["label"] + "\t"
		outstr+="\n"

	return(g, outstr)



rule gmls:
	input:
		sdaedges="PlaceOnt/SDA.{ID}.edges",
		edges = "PlaceOnt/ref.edges",
	output:
		gml="PlaceOnt/SDA.{ID}.gml",
		ref="PlaceOnt/ref.{ID}.gml",
		comps="PlaceOnt/SDA.{ID}.clusters",
		refcomps="PlaceOnt/ref.{ID}.clusters",
	run:
		realset = makeSets(input["edges"])
		Edges = pd.read_table(input["sdaedges"])
		g, outstr = makeGraph(Edges, realset)
		nx.write_gml(g, output["gml"])
		open(output["comps"], "w+").write(outstr)

		# all the ref things that doesn really matter
		Edges = pd.read_table(input["edges"])
		g, outstr = makeGraph(Edges)
		nx.write_gml(g, output["ref"])
		open(output["refcomps"], "w+").write(outstr)



rule setOverlaps:
	input:
		comps="PlaceOnt/SDA.{ID}.edges",
		refcomps="PlaceOnt/ref.edges",
	output:
		agree="PlaceOnt/agree.{ID}.txt",
	run:
		ref = makeSets(input["refcomps"])	
		sda = makeSets(input["comps"])	
		counter = 0
		for pair in sda:
			if(pair in ref):
				counter += 1
		
		per = counter*100.0/(len(sda))
		rtn = "PerCorrect:{:.02f}\tTotalCorrect:{}\tTotal:{}\n".format( per, counter/2, len(sda)/2 )
		print(rtn)
		open(output["agree"], "w+").write(rtn)



# samtools view -c PlaceOnt/OntToSDA.nof.bam ; samtools view PlaceOnt/OntToSDA.nof.bam | cut -f 3 | sort | uniq | wc -l; samtools view PlaceOnt/OntToSDA.nof.bam | cut -f 1 | sort | uniq | wc -l
#
rule final:
	input:
		sdabam = "PlaceOnt/SDAtoSDA.bam",
		gmls=expand("PlaceOnt/SDA.{ID}.gml", ID=IDS),
		agree=expand("PlaceOnt/agree.{ID}.txt", ID=IDS),
		acc=expand("PlaceOnt/AlnCorrect.{ID}.txt", ID=IDS),
		#sdaToOnt = "PlaceOnt/SDAtoOnt.bam",
		#bed = "PlaceOnt/intersect.bed",
		#bam="PlaceOnt/SDAtoRef.bam",
		#tablesbed=rules.makeTBL.output,
	output:
		final=temp("final"),
	shell:
		"touch {output}"


