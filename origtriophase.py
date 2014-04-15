#phases trios for mendelian inheritance of variants

import sys
import string
from gemini import GeminiQuery
from gemini import gemini_subjects as subjects

def phase_genotypes(database):

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
				if kid_gt_type == 0 and dad_gt_type == 1 and mom_gt_type == 1 or kid_gt_type == 3 and dad_gt_type == 1 and mom_gt_type == 1:
					mendelian = "mendelian"
					phasable = "unphasable"
					inheritance = "unknown"
					origin = "inherited from both parents"
				elif kid_gt_type == 0 and dad_gt_type == 0 and mom_gt_type == 0 or kid_gt_type == 1 and dad_gt_type == 1 and mom_gt_type == 1 or kid_gt_type == 3 and dad_gt_type == 3 and mom_gt_type == 3:
					mendelian = "mendelian"
					phasable = "unphasable"
					inheritance = "unknown"
					origin = "inherited from both parents"
				elif kid_gt_type == 2 or dad_gt_type == 2 or mom_gt_type == 2:
					mendelian = "missing allele - unknown"
					phasable = "missing allele - unknown"
					inheritance = "missing allele - unknown"
					origin = "missing allele - unknown"
				elif kid_gt_type == 1 and dad_gt_type == 0 and mom_gt_type == 0 or kid_gt_type == 1 and dad_gt_type == 3 and mom_gt_type == 3:
					mendelian = "non-mendelian"
					phasable = "unphasable"
					inheritance = "unknown"
					origin = "de novo or erroneous data"
				elif kid_gt_type == 0 and dad_gt_type == 3 and mom_gt_type == 3 or kid_gt_type == 3 and dad_gt_type == 0 and mom_gt_type == 0:
					mendelian = "non-mendelian"
					phasable = "unphasable"
					inheritance = "unknown"
					origin = "extremely rare de novo or bad data"
				elif kid_gt_type == 0 and dad_gt_type == 0 and mom_gt_type == 1 or kid_gt_type == 0 and dad_gt_type == 1 and mom_gt_type == 0 or kid_gt_type == 3 and dad_gt_type == 3 and mom_gt_type == 1 or kid_gt_type == 3 and dad_gt_type == 1 and mom_gt_type == 3:
					mendelian = "mendelian"
					phasable = "unphasable"
					inheritance = "unknown"
					origin = "inherited from both parents or unlikely de novo"
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
					inheritance = "homozygous reference from dad_gt_type, heterozygous allele from mom_gt_type"
					origin = "inherited from both parents or unlikely de novo"
					ct=0
					if m5[0] == m6[ct]:
						phasedata = m5[0]+" from dad "+m6[ct+1]+" from mom"
					else:
						phasedata = m5[0]+" from dad "+m6[ct]+" from mom"
				elif kid_gt_type == 1 and dad_gt_type == 1 and mom_gt_type == 3:
					mendelian = "mendelian"
					phasable = "phasable"
					inheritance = "heterozygous allele from dad_gt_type, homozygous alternate from mom_gt_type"
					origin = "inherited from both parents or unlikely de novo"
					ct=0
					if m5[ct] == m6[0]:
						phasedata = m5[ct+1]+" from dad "+m6[0]+" from mom"
					else:
						phasedata = m5[ct]+" from dad "+m6[0]+" from mom"
				elif kid_gt_type == 1 and dad_gt_type == 1 and mom_gt_type == 0:
					mendelian = "mendelian"
					phasable = "phasable"
					inheritance = "heterozygous allele from dad_gt_type, homozygous reference from mom_gt_type"
					origin = "inherited from both parents or unlikely de novo"
					ct=0
					if m5[ct] == m6[0]:
						phasedata = m5[ct+1]+" from dad "+m6[0]+" from mom"
					else:
						phasedata = m5[ct]+" from dad "+m6[0]+" from mom"
				elif kid_gt_type == 0 and dad_gt_type == 1 and mom_gt_type == 3 or kid_gt_type == 0 and dad_gt_type == 3 and mom_gt_type == 1 or kid_gt_type == 3 and dad_gt_type == 0 and mom_gt_type == 1 or kid_gt_type == 3 and dad_gt_type == 1 and mom_gt_type == 0:
					mendelian = "non-mendelian"
					phasable = "unphasable"
					inheritance = "unknown"
					origin = "de novo or erroneous data"
				elif kid_gt_type == 3 and dad_gt_type == 3 and mom_gt_type == 0 or kid_gt_type == 0 and dad_gt_type == 0 and mom_gt_type == 3:
					mendelian = "non-mendelian"
					phasable = "unphasable"
					inheritance = "unknown"
					origin = "de novo or erroneous data"
				elif kid_gt_type == 1 and dad_gt_type == 3 and mom_gt_type == 0:
					mendelian = "mendelian"
					phasable = "phasable"
					inheritance = "homozygous alternate from dad_gt_type, homozygous reference from mom_gt_type"
					origin = "inherited from both parents or unlikely de novo"
					phasedata = m5[0]+" from dad "+m6[0]+" from mom"
				elif kid_gt_type == 1 and dad_gt_type == 0 and mom_gt_type == 3:
					mendelian = "mendelian"
					phasable = "phasable"
					inheritance = "homozygous reference from dad_gt_type, homozygous alternate from mom_gt_type"
					origin = "inherited from both parents or unlikely de novo"
					phasedata = m5[0]+" from dad "+m6[0]+" from mom"
				elif kid_gt_type == 0 and dad_gt_type == 3 and mom_gt_type == 0 or kid_gt_type == 3 and dad_gt_type == 0 and mom_gt_type == 3:
					mendelian = "non-mendelian"
					phasable = "unphasable"
					inheritance = "unknown"
					origin = "de novo or erroneous data"
				try: 
					print chrom + "	" + start + "	" + end + "	" + family.family_id + "	" + child.name + "	" + ref + "	" + alt + "	" + gene + "	" + impact + "	" + kid_gt_ref_depths + "	" + dad_gt_ref_depths + "	" + mom_gt_ref_depths + "	" + kid_gt_alt_depths + "	" + dad_gt_alt_depths + "	" + mom_gt_alt_depths + "	" + kid_gt + "	" + dad_gt + "	" + mom_gt + "	" + str(kid_gt_type) + "	" + str(dad_gt_type) + "	" + str(mom_gt_type) + "	" + mendelian + "	" + phasable + "	" + inheritance + "	" + origin + "	" + phasedata
				except TypeError:	
					print chrom + "	" + start + "	" + end + "	" + family.family_id + "	" + child.name + "	" + ref + "	" + alt + "	" + gene + "	" + impact + "	" + kid_gt_ref_depths + "	" + dad_gt_ref_depths + "	" + mom_gt_ref_depths + "	" + kid_gt_alt_depths + "	" + dad_gt_alt_depths + "	" + mom_gt_alt_depths + "	" + kid_gt + "	" + dad_gt + "	" + mom_gt + "	" + str(kid_gt_type) + "	" + str(dad_gt_type) + "	" + str(mom_gt_type) + "	" + mendelian + "	" + phasable + "	" + inheritance + "	" + origin

phase_genotypes(sys.argv[1])