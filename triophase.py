# triophase.py by Jim Havrilla @ quinlanlab @ UVA
# phases trios for mendelian inheritance of variants

import sys
import string
import argparse
from gemini import GeminiQuery
from gemini import gemini_subjects as subjects

def phase_genotypes():

	#defines input arguments
	parser = argparse.ArgumentParser(description='Phases variants')
	parser.add_argument('-f','--input_file', default='', help='The input file; should be a SQLite .db that gemini can read')
	parser.add_argument('-p','--min_total_parent_depth', default = '0', help='The minimum total read depth for parental alleles for variants to be considered')
	parser.add_argument('-c','--min_total_child_depth', default = '0', help='The minimum total read depth for child alleles for variants to be considered')
	parser.add_argument('-m','--min_alt_parent_depth', default = '0', help='The minimum alternate read depth for parental alleles for variants to be considered')
	parser.add_argument('-a','--min_alt_child_depth', default = '0', help='The minimum alternate read depth for child alleles for variants to be considered')

	#checks minimum number of arguments
	if len(sys.argv)<1:
		parser.print_help()
		sys.exit("Where is the input file?")
    
    #parses arguments
	args = parser.parse_args()
	database=args.input_file
	mtpd=int(args.min_total_parent_depth)
	mtcd=int(args.min_total_child_depth)
	mapd=int(args.min_alt_parent_depth)
	macd=int(args.min_alt_child_depth)
	if database == '':
		sys.exit("You must supply an input file")

	gq=GeminiQuery(database)
	families = subjects.get_families(database)
	gq.run("select chrom, start, end, ref, alt, gene, impact, gts, gt_types, gt_ref_depths, gt_alt_depths from variants")
	s2i=gq.sample_to_idx
	for row in gq:
		mendelian = ""
		phasable = ""
		inheritance = ""
		origin = ""
		phasedata = ""
		chrom = str(row['chrom'])
		start = str(row['start'])
		end = str(row['end'])
		ref = str(row['ref'])
		alt = str(row['alt'])
		gene = str(row['gene'])
		impact = str(row['impact'])
		for family in families:
			dad_idx = s2i[family.father_name]
			mom_idx = s2i[family.mother_name]
			dad_gt = str(row['gts'][dad_idx])
			mom_gt = str(row['gts'][mom_idx])
			dad_gt_type = row['gt_types'][dad_idx]
			mom_gt_type = row['gt_types'][mom_idx]
			dad_gt_ref_depths = str(row['gt_ref_depths'][dad_idx])
			mom_gt_ref_depths = str(row['gt_ref_depths'][mom_idx])
			dad_gt_alt_depths = str(row['gt_alt_depths'][dad_idx])
			mom_gt_alt_depths = str(row['gt_alt_depths'][mom_idx])
			#m5=re.search('((?:\w*|\.*))/((?:\w*|\.*))',dad_gt)
			m5=string.split(dad_gt,"/")
			#m6=re.search('((?:\w*|\.*))/((?:\w*|\.*))',mom_gt)
			m6=string.split(mom_gt,"/")
			for child in family.children:
				kid_idx = s2i[str(child.name)]
				kid_gt = str(row['gts'][kid_idx])
				kid_gt_type = row['gt_types'][kid_idx]
				kid_gt_ref_depths = str(row['gt_ref_depths'][kid_idx])
				kid_gt_alt_depths = str(row['gt_alt_depths'][kid_idx])
				#code for removing variants that do not meet parameters
				if int(dad_gt_ref_depths)+int(dad_gt_alt_depths)<mtpd or int(mom_gt_ref_depths)+int(mom_gt_alt_depths)<mtpd \
				or int(kid_gt_ref_depths)+int(kid_gt_alt_depths)<mtcd \
				or int(dad_gt_alt_depths)<mapd or int(mom_gt_alt_depths)<mapd \
				or int(kid_gt_alt_depths)<macd:
					continue
				if kid_gt_type == 2 or dad_gt_type == 2 or mom_gt_type == 2:
					continue
				elif kid_gt_type == 1 and dad_gt_type == 3 and mom_gt_type == 1:
					mendelian = "mendelian"
					phasable = "phasable"
					inheritance = "homozygous alternate from dad, heterozygous allele from mom"
					origin = "inherited from both parents or unlikely de novo"
					ct=0
					if m5[0] == m6[ct]:
						phasedata = m5[0]+" from dad "+m6[ct+1]+" from mom"
					else:
						phasedata = m5[0]+" from dad "+m6[ct]+" from mom"
				elif kid_gt_type == 1 and dad_gt_type == 0 and mom_gt_type == 1:
					mendelian = "mendelian"
					phasable = "phasable"
					inheritance = "homozygous reference from dad, heterozygous allele from mom"
					origin = "inherited from both parents or unlikely de novo"
					ct=0
					if m5[0] == m6[ct]:
						phasedata = m5[0]+" from dad "+m6[ct+1]+" from mom"
					else:
						phasedata = m5[0]+" from dad "+m6[ct]+" from mom"
				elif kid_gt_type == 1 and dad_gt_type == 1 and mom_gt_type == 3:
					mendelian = "mendelian"
					phasable = "phasable"
					inheritance = "heterozygous allele from dad, homozygous alternate from mom"
					origin = "inherited from both parents or unlikely de novo"
					ct=0
					if m5[ct] == m6[0]:
						phasedata = m5[ct+1]+" from dad "+m6[0]+" from mom"
					else:
						phasedata = m5[ct]+" from dad "+m6[0]+" from mom"
				elif kid_gt_type == 1 and dad_gt_type == 1 and mom_gt_type == 0:
					mendelian = "mendelian"
					phasable = "phasable"
					inheritance = "heterozygous allele from dad, homozygous reference from mom"
					origin = "inherited from both parents or unlikely de novo"
					ct=0
					if m5[ct] == m6[0]:
						phasedata = m5[ct+1]+" from dad "+m6[0]+" from mom"
					else:
						phasedata = m5[ct]+" from dad "+m6[0]+" from mom"
				elif kid_gt_type == 1 and dad_gt_type == 3 and mom_gt_type == 0:
					mendelian = "mendelian"
					phasable = "phasable"
					inheritance = "homozygous alternate from dad, homozygous reference from mom"
					origin = "inherited from both parents or unlikely de novo"
					phasedata = m5[0]+" from dad "+m6[0]+" from mom"
				elif kid_gt_type == 1 and dad_gt_type == 0 and mom_gt_type == 3:
					mendelian = "mendelian"
					phasable = "phasable"
					inheritance = "homozygous reference from dad, homozygous alternate from mom"
					origin = "inherited from both parents or unlikely de novo"
					phasedata = m5[0]+" from dad "+m6[0]+" from mom"
				else:
					continue

				print chrom + "	" + start + "	" + end + "	" + family.family_id + "	" + child.name + "	" + ref + "	" + alt + "	" + gene + "	" + impact + "	" + kid_gt_ref_depths + "	" + dad_gt_ref_depths + "	" + mom_gt_ref_depths + "	" + kid_gt_alt_depths + "	" + dad_gt_alt_depths + "	" + mom_gt_alt_depths + "	" + kid_gt + "	" + dad_gt + "	" + mom_gt + "	" + str(kid_gt_type) + "	" + str(dad_gt_type) + "	" + str(mom_gt_type) + "	" + mendelian + "	" + phasable + "	" + inheritance + "	" + origin + "	" + phasedata

phase_genotypes()