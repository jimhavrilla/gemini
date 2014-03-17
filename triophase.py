import sys
import re
import subprocess
import fileinput
import string
from gemini import GeminiQuery

def phase_genotypes(kid_idx, mom_idx, dad_idx, gq):

	gq.run("select chrom, start, end, ref, alt, gts, gt_types from variants")
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
		kid_gt = str(row['gts'][kid_idx])
		dad_gt = str(row['gts'][dad_idx])
		mom_gt = str(row['gts'][mom_idx])
		kid_gt_type = str(row['gt_types'][kid_idx])
		dad_gt_type = str(row['gt_types'][dad_idx])
		mom_gt_type = str(row['gt_types'][mom_idx])

		#m=re.search('(chr(?:\S*\s*){5})((?:\w*|\.*)/(?:\w*|\.*),{0,1})\s*((?:\w*|\.*)/(?:\w*|\.*),{0,1})\s*((?:\w*|\.*)/(?:\w*|\.*),{0,1})\s*(\d,{0,1})(\d,{0,1})(\d,{0,1})',str(row))
		#newlist.append([m.group(1).strip(","),m.group(2).strip(","),m.group(3).strip(","),m.group(4).strip(","),m.group(5).strip(","),m.group(6).strip(","),m.group(7).strip(",")])
		#print chrom + "	" + start + "	" + end + "	" + ref + "	" + alt + "	" + kid_gt + "	" + dad_gt + "	" + mom_gt + "	" + kid_gt_type + "	" + dad_gt_type + "	" + mom_gt_type

		m4=re.search('((?:\w*|\.*))/((?:\w*|\.*))',kid_gt)
		m5=re.search('((?:\w*|\.*))/((?:\w*|\.*))',dad_gt)
		m6=re.search('((?:\w*|\.*))/((?:\w*|\.*))',mom_gt)

		# families = get_families()

		# for family in families:
		# 	kid_ = family.get_kid_idx()
		# 	dad_ = family.get_dad_idx()
		# 	mom_ = family.get_mom_idx()



		if kid_gt_type == 2 or dad_gt_type == 2 or mom_gt_type == 2:
			mendelian = "missing allele - unknown"
			phasable = "missing allele - unknown"
			inheritance = "missing allele - unknown"
			origin = "missing allele - unknown"
		elif kid_gt_type == 0 and dad_gt_type == 0 and mom_gt_type == 0 or kid_gt_type == 1 and dad_gt_type == 1 and mom_gt_type == 1 or kid_gt_type == 3 and dad_gt_type == 3 and mom_gt_type == 3:
			mendelian = "mendelian"
			phasable = "unphasable"
			inheritance = "unknown"
			origin = "inherited from both parents"
		elif kid_gt_type == 0 and dad_gt_type == 1 and mom_gt_type == 1 or kid_gt_type == 3 and dad_gt_type == 1 and mom_gt_type == 1:
			mendelian = "mendelian"
			phasable = "unphasable"
			inheritance = "unknown"
			origin = "inherited from both parents"
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
			inheritance = "homozygous alternate from dad_gt_type, heterozygous allele from mom_gt_type"
			origin = "inherited from both parents or unlikely de novo"
			ct = 1
			if m5.group(1) == m6.group(ct):
				phasedata = m5.group(1)+"from dad_gt_type"+m6.group(ct+1)+"from mom_gt_type"
			else:
				phasedata = m5.group(1)+"from dad_gt_type"+m6.group(ct)+"from mom_gt_type"
		elif kid_gt_type == 1 and dad_gt_type == 0 and mom_gt_type == 1:
			mendelian = "mendelian"
			phasable = "phasable"
			inheritance = "homozygous reference from dad_gt_type, heterozygous allele from mom_gt_type"
			origin = "inherited from both parents or unlikely de novo"
			ct = 1
			if m5.group(1) == m6.group(ct):
				phasedata = m5.group(1)+"from dad_gt_type"+m6.group(ct+1)+"from mom_gt_type"
			else:
				phasedata = m5.group(1)+"from dad_gt_type"+m6.group(ct)+"from mom_gt_type"
		elif kid_gt_type == 1 and dad_gt_type == 1 and mom_gt_type == 3:
			mendelian = "mendelian"
			phasable = "phasable"
			inheritance = "heterozygous allele from dad_gt_type, homozygous alternate from mom_gt_type"
			origin = "inherited from both parents or unlikely de novo"
			ct = 1
			if m5.group(ct) == m6.group(1):
				phasedata = m5.group(ct+1)+"from dad_gt_type"+m6.group(1)+"from mom_gt_type"
			else:
				phasedata = m5.group(ct)+"from dad_gt_type"+m6.group(1)+"from mom_gt_type"
		elif kid_gt_type == 1 and dad_gt_type == 1 and mom_gt_type == 0:
			mendelian = "mendelian"
			phasable = "phasable"
			inheritance = "heterozygous allele from dad_gt_type, homozygous reference from mom_gt_type"
			origin = "inherited from both parents or unlikely de novo"
			ct = 1
			if m5.group(ct) == m6.group(1):
				phasedata = m5.group(ct+1)+"from dad_gt_type"+m6.group(1)+"from mom_gt_type"
			else:
				phasedata = m5.group(ct)+"from dad_gt_type"+m6.group(1)+"from mom_gt_type"
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
			phasedata = m5.group(1)+"from dad_gt_type"+m6.group(1)+"from mom_gt_type"
		elif kid_gt_type == 1 and dad_gt_type == 0 and mom_gt_type == 3: #line 10 in excel sheet
			mendelian = "mendelian"
			phasable = "phasable"
			inheritance = "homozygous reference from dad_gt_type, homozygous alternate from mom_gt_type"
			origin = "inherited from both parents or unlikely de novo"
			phasedata = m5.group(1)+"from dad_gt_type"+m6.group(1)+"from mom_gt_type"
		elif kid_gt_type == 0 and dad_gt_type == 3 and mom_gt_type == 0 or kid_gt_type == 3 and dad_gt_type == 0 and mom_gt_type == 3:
			mendelian = "non-mendelian"
			phasable = "unphasable"
			inheritance = "unknown"
			origin = "de novo or erroneous data"
		try: 
			print chrom + "	" + start + "	" + end + "	" + ref + "	" + alt + "	" + kid_gt + "	" + dad_gt + "	" + mom_gt + "	" + kid_gt_type + "	" + dad_gt_type + "	" + mom_gt_type + "	" + mendelian + "	" + phasable + "	" + inheritance + "	" + origin + "	" + phasedata
		except TypeError:	
			print chrom + "	" + start + "	" + end + "	" + ref + "	" + alt + "	" + kid_gt + "	" + dad_gt + "	" + mom_gt + "	" + kid_gt_type + "	" + dad_gt_type + "	" + mom_gt_type + "	" + mendelian + "	" + phasable + "	" + inheritance + "	" + origin

gq=GeminiQuery(sys.argv[1])
s2i=gq.sample_to_idx
kid_idx = s2i['NA12878']; dad_idx = s2i['NA12891']; mom_idx = s2i['NA12892']
phase_genotypes(kid_idx,dad_idx,mom_idx,gq)