# This is a series of commands that allow creating mock files 
# from just the aggregate bcf/vcf file needed to run SiteQC pipeline

# 1. Get the list of participants from bcf/vcf file
bcftools view -h gel_mainProgramme_aggv2_chr9_67304800_67564048.vcf.gz > header.txt
tail -n 1 header.txt > line.txt
# Here manually remove all the first non-participant colnames from the header line.
tr -s '\t' '\n' < line.txt > samples.txt

# 2. Create triofile from sample list with pre-prepared header and phenotype and family values that work

awk 'BEGIN{print "trio_id\tcohort_id\tfamily_id\tinterpretation_version\tsampleId\tparticipantId\ttrio\tpedigreeId\tmotherId\tfatherId\tisProband\taffectionStatus\tsex\tpersonKaryotypicSex\tyearOfBirth\tlifeStatus\tadoptedStatus\tconsanguineousParents\tparticipantQCState\tvcf_id\tvcf_size\tgenome\tdenovo_filename"; i=3; trio[1]="Father; trio[2]="Mother"; trio[3]="Offspring"; sex[1]="MALE"; sex[2]="FEMALE"; sex[3]="FEMALE"} {i++; k=int(i/3); g=int(i%3); print "trio_" k "\tcohort1\tfamily" k "\txxxx\t" $1 "\tpid" i-2 "\t" trio[g+1] "\txxxxx\t" i*2 "\t" i*2+1 "\tFALSE\tUNAFFECTED\t" sex[g+1] "\tUNKNOWN\txxxx\tALIVE\tnotAdopted\tunknown\tpassedMedicalReviewReadyForInterpretation\t1\t3\tGRCh38\txxxx"}' samples.txt > mock_triofile_78k.tsv

