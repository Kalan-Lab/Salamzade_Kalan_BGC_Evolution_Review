import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter

def parseCDSCoord(str_gbk_loc):
	"""
	Description:
	Function to parse a string from a GenBank feature's coordinates.
	********************************************************************************************************************
	Parameters:
	- str_gbk_loc: The string value of the feature coordinate from a GenBank file being parsed by Biopython.
	********************************************************************************************************************
	"""
	try:
		start = None
		end = None
		direction = None
		all_coords = []
		is_multi_part = False
		if not 'join' in str(str_gbk_loc) and not 'order' in str(str_gbk_loc):
			start = min([int(x.strip('>').strip('<')) for x in
						 str(str_gbk_loc)[1:].split(']')[0].split(':')]) + 1
			end = max([int(x.strip('>').strip('<')) for x in
					   str(str_gbk_loc)[1:].split(']')[0].split(':')])
			direction = str(str_gbk_loc).split('(')[1].split(')')[0]
			all_coords.append([start, end, direction])
		elif 'order' in str(str_gbk_loc):
			is_multi_part = True
			all_starts = []
			all_ends = []
			all_directions = []
			for exon_coord in str(str_gbk_loc)[6:-1].split(', '):
				ec_start = min(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
				ec_end = max(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
				ec_direction = exon_coord.split('(')[1].split(')')[0]
				all_starts.append(ec_start)
				all_ends.append(ec_end)
				all_directions.append(ec_direction)
				all_coords.append([ec_start, ec_end, ec_direction])
			assert (len(set(all_directions)) == 1)
			start = min(all_starts)
			end = max(all_ends)
			direction = all_directions[0]
		else:
			is_multi_part = True
			all_starts = []
			all_ends = []
			all_directions = []
			for exon_coord in str(str_gbk_loc)[5:-1].split(', '):
				ec_start = min(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
				ec_end = max(
					[int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
				ec_direction = exon_coord.split('(')[1].split(')')[0]
				all_starts.append(ec_start)
				all_ends.append(ec_end)
				all_directions.append(ec_direction)
				all_coords.append([ec_start, ec_end, ec_direction])
			assert (len(set(all_directions)) == 1)
			start = min(all_starts)
			end = max(all_ends)
			direction = all_directions[0]
		return(all_coords, start, end, direction, is_multi_part)
	except Exception as e:
		raise RuntimeError(traceback.format_exc())

def parseAntiSMASH(bgc_genbank, comprehensive_parsing=True, flank_size=2000):
	""" Function to parse BGC Genbank produced by antiSMASH """
	bgc_info = []
	domains = []
	core_positions = set([])
	full_sequence = ""
	with open(bgc_genbank) as ogbk:
		domain_feature_types = ['PFAM_domain', 'CDS_motif', 'aSDomain']
		for rec in SeqIO.parse(ogbk, 'genbank'):
			full_sequence = str(rec.seq)
			for feature in rec.features:
				if comprehensive_parsing and feature.type in domain_feature_types:
					all_coords, start, end, direction, is_multi_part = parseCDSCoord(str(feature.location))
					aSDomain = "NA"
					description = "NA"
					try:
						aSDomain = feature.qualifiers.get('aSDomain')[0]
					except:
						pass
					try:
						description = feature.qualifiers.get('description')[0]
					except:
						pass
					domains.append({'start': start, 'end': end, 'type': feature.type, 'aSDomain': aSDomain,
									'description': description, 'is_multi_part': is_multi_part})
				elif feature.type == 'protocluster':
					detection_rule = feature.qualifiers.get('detection_rule')[0]
					try:
						product = feature.qualifiers.get('product')[0]
					except:
						product = "NA"
					contig_edge = feature.qualifiers.get('contig_edge')[0]
					bgc_info.append(
						{'detection_rule': detection_rule, 'product': product, 'contig_edge': contig_edge,
						 'full_sequence': str(rec.seq)})
				elif feature.type == 'proto_core':
					if not 'join' in str(feature.location):
						core_start = min([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
						core_end = max([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
						core_positions = core_positions.union(set(range(core_start + 1, core_end + 1)))
					else:
						core_starts = []
						core_ends = []
						for exon_coord in str(feature.location)[5:-1].split(', '):
							start = min([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
							end = max([int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
							core_starts.append(start); core_ends.append(end)
						core_positions = core_positions.union(set(range(min(core_starts)+1, max(core_ends)+1)))

	if len(bgc_info) == 0:
		bgc_info = [{'detection_rule': 'NA', 'product': 'NA', 'contig_edge': 'NA', 'full_sequence': full_sequence}]

	genes = {}
	core_genes = set([])
	gene_order = {}
	with open(bgc_genbank) as ogbk:
		for rec in SeqIO.parse(ogbk, 'genbank'):
			for feature in rec.features:
				if feature.type == "CDS":
					lt = feature.qualifiers.get('locus_tag')[0]
					all_coords, start, end, direction, is_multi_part = parseCDSCoord(str(feature.location))

					try:
						product = feature.qualifiers.get('product')[0]
					except:
						product = "hypothetical protein"

					rule_based_bgc_cds = False
					try:
						if 'rule-based-clusters' in feature.qualifiers.get('gene_functions')[0]:
							rule_based_bgc_cds = True
					except:
						pass

					grange = set(range(start, end + 1))
					core_overlap = False
					if len(grange.intersection(core_positions)) > 0 and rule_based_bgc_cds:
						core_overlap = True
						core_genes.add(lt)

					gene_order[lt] = start

					prot_seq, nucl_seq, nucl_seq_with_flanks, relative_start, relative_end, gene_domains = [None] * 6
					if comprehensive_parsing:
						prot_seq = feature.qualifiers.get('translation')[0]
						gene_domains = []
						for d in domains:
							drange = set(range(d['start'], d['end'] + 1))
							if len(drange.intersection(grange)) > 0:
								gene_domains.append(d)

						flank_start = start - flank_size
						flank_end = end + flank_size

						if flank_start < 1: flank_start = 1

						if flank_end >= len(full_sequence): flank_end = None
						if end >= len(full_sequence): end = len(full_sequence)

						nucl_seq = ''
						for sc, ec, dc in sorted(all_coords, key=itemgetter(0), reverse=False):
							if ec >= len(full_sequence):
								nucl_seq += full_sequence[sc - 1:]
							else:
								nucl_seq += full_sequence[sc - 1:ec]

						if flank_end:
							nucl_seq_with_flanks = full_sequence[flank_start - 1:flank_end]
						else:
							nucl_seq_with_flanks = full_sequence[flank_start - 1:]

						gene_length = end - start

						relative_start = nucl_seq_with_flanks.find(nucl_seq)
						relative_end = relative_start + gene_length

						if direction == '-':
							nucl_seq = str(Seq(nucl_seq).reverse_complement())
							nucl_seq_with_flanks = str(Seq(nucl_seq_with_flanks).reverse_complement())
							relative_start = nucl_seq_with_flanks.find(nucl_seq)
							relative_end = relative_start + gene_length

					genes[lt] = {'start': start, 'end': end, 'direction': direction,
								 'product': product, 'prot_seq': prot_seq, 'nucl_seq': nucl_seq,
								 'nucl_seq_with_flanks': nucl_seq_with_flanks, 'gene_domains': gene_domains,
								 'core_overlap': core_overlap, 'relative_start': relative_start,
								 'relative_end': relative_end, 'is_multi_part': is_multi_part}

	number_of_core_gene_groups = 0
	tmp = []
	for lt in sorted(gene_order.items(), key=itemgetter(1), reverse=True):
		if lt[0] in core_genes:
			tmp.append(lt[0])
		elif len(tmp) > 0:
			number_of_core_gene_groups += 1
			tmp = []
	if len(tmp) > 0:
		number_of_core_gene_groups += 1

	for i, pc in enumerate(bgc_info):
		bgc_info[i]['count_core_gene_groups'] = number_of_core_gene_groups

	return(genes)

antismash_results = os.path.abspath('lsaBGC_Results/Recreated_GenBanks_with_LocusTags/') + '/'

gene_to_core_status = {}
for g in os.listdir(antismash_results):
	gca_dir = antismash_results + g + '/'
	for f in os.listdir(gca_dir):
		if f.endswith('.gbk') and 'region' in f:
			gbk_file = gca_dir + f 
			genes = parseAntiSMASH(gbk_file)
			for gene in genes:
				gene_to_core_status['.'.join(f.split('.')[:-1]) + '|' + gene] = genes[gene]['core_overlap']
				#print(g + '\t' + '.'.join(f.split('.')[:-1]) + '\t' + f + '\t' + gene + '\t' + str(genes[gene]['core_overlap']))

with open('lsaBGC_Results/Scratch_Workspace/zol_results.tsv') as ozo:
	for i, line in enumerate(ozo):
		if i == 0: continue
		line = line.strip('\n')
		ls = line.split('\t')
		gcf, og, single_copy, tot_prop, comp_prop, med_len = ls[:6]
		tajd = ls[8]
		tot_prop = float(tot_prop)
		lts = ls[-2].split('; ')
		if single_copy == 'True' and tot_prop >= 0.80 and len(lts) >= 4:
			core_count = 0
			tot_count = 0
			for lt in lts:
				tot_count += 1
				if gene_to_core_status[lt]:
					core_count += 1

			protocore = 'False'
			if (core_count/float(tot_count)) >= 0.5:
				protocore = 'True'
			if tajd != 'NA' and not tajd.startswith("< 3"):	
				print(gcf + '\t' + og + '\t' + protocore + '\t' + tajd + '\t' + med_len + '\t' + str(len(lts)) + '\t' + str(tot_prop))
		
