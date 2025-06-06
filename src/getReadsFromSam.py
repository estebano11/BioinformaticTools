import pysam
import argparse
import os
import sys
from multiprocessing import Pool, cpu_count

# --------------------*********---------------------
# Script's description
parser = argparse.ArgumentParser(description = "Take all the reads mapped from a MAG")
parser.add_argument('-f', '--mag',dest='mag', help="mag file",required=True)
parser.add_argument('-e', '--extension',dest='ext', default='fasta', help="MAG extension")
parser.add_argument('-S', '--sam', dest='sam', nargs='+', help="Sam file, if is more than one a comma separated list or just with *sam",required=True)
parser.add_argument('-q','--quality', dest='quality', default='any',help="caz-coverage",required=True)
parser.add_argument('-t','--threads', dest='threads', default=1, type=int, help="caz-coverage")
parser.add_argument('-o', '--outdir', dest='out', help="Folder to store the fastq Files",required=True)
args = parser.parse_args()
# --------------------*********---------------------
print ("--------------*****------------------")
# --------------------*********---------------------
def readSamFile(sam_file):
	samfile = pysam.AlignmentFile(sam_file, 'rb')
	allReadsDict = {read.query_name: read for read in samfile}  # Store reads in a dict for fast lookup
	allReadsSet = set(allReadsDict.keys())  # Use the dictionary keys for set operations
	return samfile, allReadsSet, allReadsDict

# --------------------*********---------------------
def findReads(sam_file,out_dir,out_file,contigs,quality):
	samName = str(sam_file).replace(".sorted.bam","").split("/")[-1]
	print (f"Reading {samName} file .... ")
	samfile, allReadsSet, allReadsDict = readSamFile(sam_file)
	print (f"Finished reading bam file {samName}")
	mappedReads = pysam.AlignmentFile(f"{out_dir}/{out_file}.mappedReads.bam", "wb", template=samfile)
	
	fastq_forward = open(f"{out_dir}/{out_file}.forward.fastq", "w")
	fastq_reverse = open(f"{out_dir}/{out_file}.reverse.fastq", "w")
	fastq_single = open(f"{out_dir}/{out_file}.single.fastq", "w")
	
	pairedReads = {}
	single = set()
	readsMapped = []
	
	# Step 1: Collect reads and track pairing
	for contig in contigs:
		for read in samfile.fetch(contig):
			if read.mapping_quality >= quality or read.mapping_quality == 0:
				readName = read.query_name
				readsMapped.append(readName)
				readN = readName.split("_")[0]  # Extract base read name (without _1/_2)
				# Add to single list (will remove it when the pair is found)
				if readN not in single:
					single.add(readN)
				else:
					single.remove(readN)  # Remove if we already have the pair

				# Initialize in pairedReads if not already present
				if readN not in pairedReads:
					pairedReads[readN] = [None, None]

				# Assign read to forward (_1) or reverse (_2) pair
				if "_1" in readName:
					forward = read
					reverse = pairedReads[readN][1]
					pairedReads[readN] = [forward, reverse]
				elif "_2" in readName:
					forward = pairedReads[readN][0]
					reverse = read
					pairedReads[readN] = [forward, reverse]
				else:
					print("Your read doesn't have identification (_1 or _2)")
					sys.exit(123)

	# Step 2: Search for missing pairs in `single`
	# allReadsSet = set(pairedReads.keys())  # Convert pairedReads keys to a set
	print(f"Initial unpaired reads count: {len(single)} in {samName}")

	for readN in list(single):
		forward, reverse = pairedReads[readN]
		if forward is None:  # If forward is missing
			f = reverse.query_name.replace("_2", "_1")
			if f in allReadsSet:
				fwd_read = allReadsDict[f]
				pairedReads[read] = [fwd_read, reverse]
				readsMapped.append(fwd_read)
				single.remove(readN)

		elif reverse is None:  # If reverse is missing
			r = forward.query_name.replace("_1", "_2")
			if r in allReadsSet:
				rev_read = allReadsDict[r]
				pairedReads[read] = [forward, rev_read]
				readsMapped.append(rev_read)
				single.remove(readN)

	print(f"Remaining unpaired reads after search: {len(single)} in {samName}")
	print(f"Total mapped reads: {len(readsMapped)} in {samName}")
	# Step 3: Write paired and unpaired reads to the BAM file
	for readN, pair in pairedReads.items():
		forward, reverse = pair

		# If both reads are available, write them both
		if forward is not None and reverse is not None:
			mappedReads.write(forward)
			mappedReads.write(reverse)

			write_fastq_entry(fastq_forward, forward)
			write_fastq_entry(fastq_reverse, reverse)
		# If only one of the pair exists, write the one that exists
		elif forward is not None:
			mappedReads.write(forward)
			write_fastq_entry(fastq_single, forward)
		elif reverse is not None:
			mappedReads.write(reverse)
			write_fastq_entry(fastq_single, reverse)

	# Close the BAM output file
	mappedReads.close()
	samfile.close()
	fastq_forward.close()
	fastq_reverse.close()
	fastq_single.close()
	print(f"BAM file written in {out_dir}/{out_file}.mappedReads.bam")

# -------------*****------------------
# Helper function to write a FASTQ entry
def write_fastq_entry(fastq_file, read):
	# Extract the fields for FASTQ format
	read_name = f"@{read.query_name}"
	sequence = read.query_sequence
	plus_line = "+"
	# Convert quality scores to ASCII format (FASTQ uses Phred+33 encoding)
	qualities = "".join(chr(q + 33) for q in read.query_qualities)
	# Write the FASTQ entry
	fastq_file.write(f"{read_name}\n{sequence}\n{plus_line}\n{qualities}\n")

# -------------*****------------------
def process_sam_file(sam_file,contigs):
	out_dir = os.path.abspath(args.out)
	extension = args.ext
	quality = int(args.quality)
	out_file = magsFile.replace('.' + extension, '').split("/")[-1]
	sample = sam_file.replace(".sorted.bam", "").split("/")[-1]
	out_file = f"{out_file}_{sample}"
	out_dir = os.path.abspath(out_dir)

	findReads(sam_file, out_dir, out_file, contigs, quality)
# -------------*****------------------
threads = int(args.threads)
magsFile = args.mag
out_dir = os.path.abspath(args.out)
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
	print("Directory created successfully!")
else:
	print("Directory already exists!")
# Handle multiple SAM files with wildcard (*) or comma-separated list
if "," in args.sam:
		sam_files = args.sam.split(",")
else:
		sam_files = args.sam

if len(sam_files) == 0:
		print("No SAM files found. Exiting.")
		sys.exit(1)
print (f"Analizing mapped reads in {len(sam_files)} BAM files")
# print (f"{sam_files}")

# -------------*****------------------
# Primero me quedo con los nombres de los contigs
contigs = []
with open(magsFile)as f:
	for line in f:
		if line.startswith(">"):
			contig = line.strip().replace(">","")
			contigs.append(contig)
print (f'Your MAG has {len(contigs)} contigs')
# -------------*****------------------
# Run in parallel if more than one thread is specified
if threads > 1:
	print(f"Running with {threads} threads")
	# Prepare tuples of arguments for starmap
	args_for_pool = [(sam_file, contigs) for sam_file in sam_files]
	with Pool(threads) as pool:
		pool.starmap(process_sam_file, args_for_pool)
else:
	print ("Running iteratively with one core")
	total_samples = len(sam_files)
	x = 1
	for sam_file in sam_files:
		print (f"Running sample {x} out of {total_samples}")
		process_sam_file(sam_file, contigs)
		x += 1