
#!/bin/bash
#
# @ Author: Daniel Rhodes
# @ Create Time: 2020-09-15 11:16:12
# @ Description: Site QC metrics for the aggV2
#


#This is the version being created as a cleaned up version of aggV2 run. Actual original code is found in sample_metrics.sh, which
#serves as an indicator of how the code was actually run. This is an effort to tidy up the code an optimise.
#TODO
#Optimise median calcs by removing sed command and incorporating into awk

#Site QC functions

#Define trios
triodata_define(){
    ###############
    ##  Purpose: Create a .fam file of confirmed trios in our samples
    ##  Input: List of trios, samples included in agg
    ##  Output: .fam and .keep files
    ###############
    module load $RLoad
    Rscript trio_define.R $triodata $aggregateSamples $triodata.fam $triodata.keep
}

completeSites(){
    ###############
    ##  Purpose: Make sure the number of samples is listed in resources
    ##  Input: BCF
    ##  Output: Txt file with N samples
    ###############
    module load $bcftoolsLoad
    export infile=`sed -n "1p" $bedfile`
    #All this will do is create a file containing the number of samples to be used
    #by the annotation script
    if [ ! -f "${resources}/N_samples" ]; then
        bcftools query -l ${input}${infile} | wc -l > ${resources}/N_samples
    fi
}

startFile(){
    ###############
    ##  Purpose: Create a backbone of IDs for other data to be joined to
    ##  Input: BCF
    ##  Output: txt file with IDs
    ###############
    module load $bcftoolsLoad
    mkdir -p ${out}startfile
    echo 'Creating backbone file'
    export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
    if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
    
    # Query chrom and pos for annotation file
    if "$sexChrom"; then
        #XY
        bcftools query -f '%CHROM %POS %REF %ALT %INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC\n' ${input}{${infile} -S ${resources}${xy} --output ${out}startfile/start_file_${i}_XY
            #XX
            bcftools query -f '%CHROM %POS %REF %ALT %INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC\n' ${input}${infile} -S ${resources}${xx} --output ${out}startfile/start_file_${i}_XX
        else
            #Autosomes
            bcftools query -f '%CHROM %POS %REF %ALT %INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC\n' ${input}${infile} --output ${out}startfile/start_file_${i}
        fi
    }
    
    missingness1(){
        ###############
        ##  Purpose: Missingness step 1, count fully missing GTs
        ##  Input: BCF
        ##  Output: Txt file with ID and count
        ###############
        module load $bcftoolsLoad
        echo 'Calculating missing sites'
        mkdir -p ${out}missing
        module load $bcftoolsLoad
        
        export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
        export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
        if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
        if "$sexChrom"; then
            #XY
            bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -i 'GT="./." & FORMAT/DP=0'  -S ${resources}${xy} \
            | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing1_${i}_XY
            #XX
            bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -i 'GT="./." & FORMAT/DP=0'  -S ${resources}${xx} \
            | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing1_${i}_XX
        else
            #autosomes
            bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -i 'GT="./." & FORMAT/DP=0' \
            | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing1_${i}
        fi
    }
    
    missingness2(){
        ###############
        ##  Purpose: Count complete GTs only
        ##  Input: BCF
        ##  Output: Txt file with ID and count
        ###############
        module load $bcftoolsLoad
        echo 'Calculating missing sites'
        mkdir -p ${out}missing
        module load $bcftoolsLoad
        
        export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
        export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
        if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
        if "$sexChrom"; then
            #XY
            bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -e 'GT~"\."' -S ${resources}${xy} \
            | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing2/missing2_${i}_XY
            #XX
            bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -e 'GT~"\."' -S ${resources}${xx} \
            | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing2_${i}_XX
        else
            #Autosomes
            bcftools query ${input}${infile} -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GT ] \n' -e 'GT~"\."' \
            | awk 'BEGIN{FS=" "} {print $1" "NF-1}' > ${out}missing/missing2_${i}
        fi
    }
    
    medianCovAll(){
        ###############
        ##  Purpose: Produce median value for depth across all GT
        ##  Input: BCF
        ##  Output: Txt file with ID and median depth
        ###############
        module load $bcftoolsLoad
        mkdir -p ${out}/medianCovAll
        echo 'Calculating median depth for all GT...'
        export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
        export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
        if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
        
        if "$sexChrom"; then
            #XY
            bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' ${input}${infile} -S ${resources}${xy} |  \
            sed s/[[:space:]]\\./\ 0/g | \
            awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
            print $1"\t"median}' > ${out}medianCovAll/medianCov_${i}
            #XX
            bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' ${input}${infile} -S ${resources}${xx}|  \
            sed s/[[:space:]]\\./\ 0/g | \
            awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
            print $1"\t"median}' > ${out}medianCovAll/medianCov_${i}
        else
            #Autosomes
            bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' ${input}${infile} |  \
            sed s/[[:space:]]\\./\ 0/g | \
            awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
            print $1"\t"median}' > ${out}medianCovAll/medianCov_${i}
        fi
    }
    
    medianCovNonMiss(){
        module load $bcftoolsLoad
        mkdir -p ${out}/medianCovNonMiss
        echo 'Calculating median depth for non missing GT...'
        export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
        export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
        if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
        
        if "$sexChrom"; then
            #XY
            bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' -e 'GT~"\."' ${input}${infile} -S ${resources}${xy} |  \
            sed s/[[:space:]]\\./\ 0/g | \
            awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
            print $1"\t"median}' > ${out}medianCovNonMiss/medianCovNonMiss_${i}
            #XX
            bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' -e 'GT~"\."' ${input}${infile} -S ${resources}${xx}|  \
            sed s/[[:space:]]\\./\ 0/g | \
            awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
            print $1"\t"median}' > ${out}medianCovNonMiss/medianCovNonMiss_${i}
        else
            #Autosomes
            bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [ %DP]\n' -e 'GT~"\."' ${input}${infile} |  \
            sed s/[[:space:]]\\./\ 0/g | \
            awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2;
            print $1"\t"median}' > ${out}medianCovNonMiss/medianCovNonMiss_${i}
        fi
    }
    
    medianGQ(){
        ###############
        ##  Purpose: Calculate median GQ
        ##  Input: BCF
        ##  Output: txt file - ID and median
        ###############
        module load $bcftoolsLoad
        echo 'Calculating median GQ for non missing GT...'
        export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
        export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
        if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
        mkdir -p ${out}/medianGQ
        
        if "$sexChrom"; then
            #XY
            bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GQ ] \n' -e 'GT~"\."' ${input}${infile} -S ${resources}${xy} |\
            sed s/[[:space:]]\\./\ 0/g | \
            awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2; if(median > 99) {median=99};
            print $1"\t"median}' > ${out}medianGQ/medianGQ_${i}_XY
            #XX
            bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GQ ] \n' -e 'GT~"\."' ${input}${infile} -S ${resources}${xx} |\
            sed s/[[:space:]]\\./\ 0/g | \
            awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2; if(median > 99) {median=99};
            print $1"\t"median}' > ${out}medianGQ/medianGQ_${i}_XX
        else
            bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%GQ ] \n' -e 'GT~"\."' ${input}${infile} |\
            sed s/[[:space:]]\\./\ 0/g | \
            awk -F '[[:space:]]+' '{ n=split($0,a)-1; asort(a);
                median=n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2; if(median > 99) {median=99};
            print $1"\t"median}' > ${out}medianGQ/medianGQ_${i}
        fi
    }
    
    ABRatioP1(){
        ###############
        ##  Purpose: AB ratio calculation - number of hets passing binomial test (reads supporting het call)
        ##  Input: BCF
        ##  Output: Txt file with Nhets that pass
        ###############
        module load $bcftoolsLoad
        echo 'AB ratio part 1'
        mkdir -p ${out}AB_hetPass
        export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
        export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
        if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
        
        if "$sexChrom"; then
            #We only calculate AB ratio for XX
            bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n'  -S ${resources}${xx} \
            -i 'GT="het" & binom(FMT/AD) > 0.01' ${input}${infile} | \
            awk '{print $1"\t"NF -1}' > ${out}AB_hetPass/hetPass_${i}_XX
        else
            bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n' \
            -i 'GT="het" & binom(FMT/AD) > 0.01' ${input}${infile} | \
            awk '{print $1"\t"NF -1}' > ${out}AB_hetPass/hetPass_${i}
        fi
    }
    
    ABRatioP2(){
        ###############
        ##  Purpose: Number of het GTs for p2 AB ratio
        ##  Input: BCF
        ##  Output: txt file with ID and N hets
        ###############
        module load $bcftoolsLoad
        echo 'AB ratio part 2'
        mkdir -p ${out}AB_hetAll
        export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
        export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
        if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
        
        if "$sexChrom"; then
            #We only calculate AB ratio for XX
            bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n' -i 'GT="het"' ${input}${infile} -S ${resources}${xx} | \
            awk '{print $1"\t"NF-1}' > ${out}AB_hetAll/hetAll_${i}_XX
        else
            bcftools query -f '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC [%AD ] \n' -i 'GT="het"' ${input}${infile}  | \
            awk '{print $1"\t"NF-1}' > ${out}AB_hetAll/hetAll_${i}
        fi
    }
    
    pull_1KGPSites(){
        ###############
        ##  Purpose: Pull sites from 1000KGP3
        ##  Input: 1000KGP3 genotypes vcf
        ##  Output: Txt file
        ###############
        mkdir -p ${out}/1KGP3_intersect/
        chr="${LSB_JOBINDEX}"
        kgpin=$(ls -d "/public_data_resources/1000-genomes/20130502_GRCh38/"*.vcf.gz | \
        grep ${chr}_ | grep -v genotypes )
        
        zcat ${kgpin} | awk '/^##/ {next} { print $1"\t"$2"\t"$4"\t"$5}'  > tmp_1kgp_${chr}.txt
    }
    
    aggregateAnnotation(){
        ###############
        ##  Purpose: Annotate and make pass/fail. If king set to T in env, print subset of cols
        ##  Input: All the outputs of step 1 metrics
        ##  Output: A single text file containing annotated variants
        ###############
        #TODO - add option for chrom X here
        module load $RLoad
        
        mkdir -p ${out}Annotation
        #Location for summary data
        mkdir -p ${out}Annotation/Summary_stats #loc hardcoded in R script
        export infile=`sed -n "${LSB_JOBINDEX}p" $bedfile`
        export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
        if [[ "$infile" == *"chrX"* ]]; then sexChrom=true; fi
        
        export R_LIBS="$(dirname `which R`)/../lib64/R/library" #Set me to mitigate library mismatch on helix
        #Remember to set path var for .libpath
        Rscript annotatePerChunk.R \
        ${i} \
        ${out}startfile/start_file_${i} \
        ${out}missing/missing1_${i} \
        ${out}missing/missing2_${i} \
        ${out}medianCovAll/medianCov_${i} \
        ${out}medianCovNonMiss/medianCovNonMiss_${i} \
        ${out}medianGQ/medianGQ_${i} \
        ${out}AB_hetAll/hetAll_${i} \
        ${out}AB_hetPass/hetPass_${i} \
        ${out}MendelErrSites/MendErr_${i}.lmendel \
        ${out}Annotation \
        ${resources}/N_samples \
        ${out}AC_counts/${i}_AC \
        tmp_1kgp_${chr}.txt \
        ignore #placeholder for hardy
    }
    
    
    ### KING WORKFLOW ###
    
    sort_compress(){
        ###############
        ##  Purpose: Sort and compress site metric data for KING step
        ##  Input: Text file containing annotated variants
        ##  Output: bgzipped and tabixed metrics information
        ###############
        
        autosomes=`grep -v chrX ${bedfile} | grep -v chrY | grep -v chrX | grep -v chrM`
        export infile=`echo $autosomes | cut -d' ' -f ${LSB_JOBINDEX}`
        #Only autosomes
        export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
        
        module load $bcftoolsLoad
        sort -k2 -n ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}.txt > \
        ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt
        bgzip -f ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt && \
        tabix -s1 -b2 -e2 ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt.gz
        rm ${out}Annotation/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt
    }
    
    regions_filter(){
        ###############
        ##  Purpose: Produce BCFs of our data filtered to sites pass sites
        ##  Input: BCF, site metrics KING sites
        ##  Output: filtered bcf
        ###############
        module load $bcftoolsLoad
        mkdir -p ${out}AnnotatedVCFs/regionsFiltered
        autosomes=`grep -v chrX ${bedfile} | grep -v chrY | grep -v chrX | grep -v chrM`
        export infile=`echo $autosomes | cut -d' ' -f ${LSB_JOBINDEX}`
        #Only autosomes
        export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
        
        bcftools view ${input}${infile} \
        -T ${out}Annotation_newtest/BCFtools_site_metrics_SUBCOLS${i}_sorted.txt.gz  \
        -Ob \
        -o ${out}AnnotatedVCFs/regionsFiltered/${i}_regionsFiltered.bcf
    }
    
    #Regions_filter and further_filtering split for speed - split bcftools runs into multiple ops
    
    further_filtering(){
        ###############
        ##  Purpose: Second stage filtering to give biallelic SNPs intersected with 1000KGP3 with MAF > 0.01
        ##  Input: Filtered BCFs, Michigan high LD pos
        ##  Output: BCF file, and txt file with positions
        ###############
        #Previously used to filter out Michigan high LD sites in this func, this now occurs
        #downstream using plink
        module load $bcftoolsLoad
        autosomes=`grep -v chrX ${bedfile} | grep -v chrY | grep -v chrX | grep -v chrM`
        export infile=`echo $autosomes | cut -d' ' -f ${LSB_JOBINDEX}`
        #Only autosomes
        export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
        
        bcftools view ${out}AnnotatedVCFs/regionsFiltered/${i}_regionsFiltered.bcf \
        -i 'INFO/OLD_MULTIALLELIC="." & INFO/OLD_CLUMPED="."' \
        -v snps  | \
        bcftools annotate \
        --set-id '%CHROM:%POS-%REF/%ALT-%INFO/OLD_CLUMPED-%INFO/OLD_MULTIALLELIC' | \
        bcftools +fill-tags -Ob \
        -o ${out}AnnotatedVCFs/regionsFiltered/MichiganLD_regionsFiltered_${i}.bcf \
        -- -t MAF
        #Produce filtered txt file
        bcftools query ${out}AnnotatedVCFs/regionsFiltered/MichiganLD_regionsFiltered_${i}.bcf \
        -i 'MAF[0]>0.01' -f '%CHROM\t%POS\t%REF\t%ALT\t%MAF\n' | \
        awk -F "\t" '{ if(($5 == "G" && $6 == "C") || ($6 == "G" && $5 == "C") || ($5 == "A" && $6 == "T") || ($6 == "A" && $5 == "T")) {next} { print $0} }' \
        > ${out}AnnotatedVCFs/MAF_filtered_1kp3intersect_${i}.txt
    }
    
    final_KING_BCF(){
        ###############
        ##  Purpose: Produce new BCF just with filtered sites
        ##  Input: Txt file of samples in agg, filtered regions bcf, intersected sites txt file
        ##  Output: Filtered compressed and indexed vcf
        ###############
        mkdir -p ${out}/KING
        autosomes=`grep -v chrX ${bedfile} | grep -v chrY | grep -v chrX | grep -v chrM`
        export infile=`echo $autosomes | cut -d' ' -f ${LSB_JOBINDEX}`
        #Only autosomes
        export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
        #Now filter down our file to just samples we want in our GRM. This removes any withdrawals that we learned of during the process of aggregation
        #Store the header
        bcftools view \
        -S ${resources}${sampleList} \
        --force-samples \
        -h ${out}AnnotatedVCFs/regionsFiltered/MichiganLD_regionsFiltered_${i}.bcf \
        > ${out}KING/${i}_filtered.vcf
        
        #Then match against all variant cols in our subsetted bcf to our maf filtered, intersected sites and only print those that are in the variant file.
        #Then append this to the stored header, SNPRelate needs vcfs so leave as is
        bcftools view \
        -H ${out}AnnotatedVCFs/regionsFiltered/MichiganLD_regionsFiltered_${i}.bcf \
        -S ${resources}${sampleList} \
        --force-samples \
        | awk -F '\t' 'NR==FNR{c[$1$2$3$4]++;next}; c[$1$2$4$5] > 0' ${out}AnnotatedVCFs/MAF_filtered_1kp3intersect_${i}.txt - >> ${out}KING/${i}_filtered.vcf
        bgzip ${out}KING/${i}_filtered.vcf
        tabix ${out}KING/${i}_filtered.vcf.gz
    }
    
    concat_KING_VCF(){
        ###############
        ##  Purpose: Concatenate compressed vcfs to per chromosome files
        ##  Input: Compressed vcfs
        ##  Output: 1 per chrom compressed vcf
        ###############
        module load $bcftoolsLoad
        i="${LSB_JOBINDEX}"
        mkdir -p ${out}perChrom_KING
        
        find ${out}KING -type f -name "chr${i}_*.vcf.gz" > tmp.files_chrom${i}.txt
        bcftools concat \
        -f tmp.files_chrom${i}.txt \
        -Oz \
        -o ${out}perChrom_KING/chrom${i}_merged_filtered.vcf.gz && \
        tabix ${out}perChrom_KING/chrom${i}_merged_filtered.vcf.gz && \
        rm tmp.files_chrom${i}.txt
    }
    
    # POSSIBLE CHECK - can read file and check for ACTG combos that should have been filtered
    makeBedAll(){
        ###############
        ##  Purpose: Make BED files for 1000KGP3 intersected vcfs
        ##  Input: per chrom vcf.gz
        ##  Output: BED files per chrom
        ###############
        module load $plink2Load
        module load $bcftoolsLoad
        echo 'Creating bed file'
        mkdir -p ${out}BEDref
        module load $plink2Load
        autosomes=`grep -v chrX ${bedfile} | grep -v chrY | grep -v chrX | grep -v chrM`
        export infile=`echo $autosomes | cut -d' ' -f ${LSB_JOBINDEX}`
        #Only autosomes
        export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
        
        mkdir -p ${out}BED
        plink2 --vcf ${out}perChrom_KING/chrom${i}_merged_filtered.vcf.gz \
        --make-bed \
        --vcf-half-call m \
        --set-missing-var-ids chr@:#-\$r/\$a-.-. \
        --new-id-max-allele-len 60 missing \
        --exclude range ${resources}${michiganList} \
        --double-id \
        --real-ref-alleles \
        --allow-extra-chr \
        --threads 30 \
        --out ${out}BED/BED_${i}
    }
    
    LD_BED(){
        ###############
        ##  Purpose: LD prune SNPs
        ##  Input: Previoulsy produced BED, BIM, FAM files
        ##  Output: BED without high LD SNPs
        ###############
        #Not considering founders in this as all of our SNPs are common
        module load $plinkLoad
        plink  \
        --keep-allele-order \
        --bfile ${out}BED/BED_${i} \
        --indep-pairwise 500kb 1 0.1 \
        --threads 30 \
        --out ${out}BED/BED_LD_${i}
        
        #Now that we have our correct list of SNPs (prune.in), filter the original
        #bed file to just these sites
        plink \
        --make-bed \
        --bfile ${out}BED/BED_${i} \
        --keep-allele-order \
        --extract ${out}/BED/BED_LD_${i}.prune.in \
        --double-id \
        --allow-extra-chr \
        --threads 30 \
        --out ${out}BED/BED_LDpruned_${i}
    }
    
    merge_autosomes(){
        ###############
        ##  Purpose: Merge autosomes to genome wide BED files
        ##  Input: per chrom bed files (LD pruned)
        ##  Output: genome wide BED files
        ###############
        module load $plinkLoad
        for i in {1..22}; do
            echo ${out}BED/BED_LDpruned_$i >> mergelist.txt
        done
        plink --merge-list mergelist.txt \
        --make-bed \
        --out ${out}BED/autosomes_LD_pruned_1kgp3Intersect
        rm mergelist.txt
    }
    
    hwe_pruning_30k_snps()){
    ###############
    ##  Purpose: produce a first pass HWE filter.
    ##  Input: LD pruned snps, bed files
    ##  Output: per pop unrelated bed files, list of hwe filtered snps
    ###############
    # Note the input file for unrels (third supplied argument to R) MUST have a column
    # 'unrelated_set' with either 1 or 0 denoting unrelated (1) or related (0)
    # Currently I have no provision to test for this
    module load $RLoad
    module load $plinkLoad
    R -e 'library(data.table);
    library(dplyr);
    args <- commandArgs(trailingOnly = T);
    setwd(args[1]);
    dat <- fread(args[2]) %>% as_tibble();
    unrels <- fread(args[3]) %>% as_tibble() %>% filter(unrelated_set == 1);
    dat <- dat %>% filter(plate_key %in% unrels$plate_key);
    for(col in c("AFR","EUR","SAS","EAS")){dat[dat[col]>0.8,c("plate_key",col)] %>% write.table(paste0(col,"pop.txt"), quote = F, row.names=F, sep = "\t")}' \
    --args \
    ${out}BED/ \
    ${resources}${snps30k} \
    ${resources}${pcs30k}
    
    #Include function check - if the pop.txt file wc -l == 1, don't include in loop,
    #this could just be added to the for loop
    #CHECK THAT THIS WILL WRITE TO THE CORRECT PLACE
    bedmain="${out}/BED/autosomes_LD_pruned_1kgp3Intersect"
    for pop in AFR EUR SAS EAS; do
        echo ${pop}
        awk '{print $1"\t"$1}' ${pop}pop.txt > ${pop}keep
        plink \
        --keep ${pop}keep \
        --keep-allele-order \
        --make-bed \
        --bfile ${bedmain} \
        --out ${pop}
        
        plink --bfile ${pop} --hardy midp --out ${pop} --nonfounders
    done
    
    
    #Combine the HWE and produce a list of pass
    R -e 'library(data.table);
    library(dplyr);
    args <- commandArgs(trailingOnly = T);
    setwd(args[1]);
    dat <- lapply(c("EUR.hwe","AFR.hwe", "SAS.hwe", "EAS.hwe"),fread);
    names(dat) <- c("EUR.hwe","AFR.hwe", "SAS.hwe", "EAS.hwe");
    dat <- dat %>% bind_rows(.id="id");
    write.table(dat, "combinedHWE.txt", row.names = F, quote = F)
    #Create set that is just SNPS that are >1e-5 in all pops
    dat %>% filter(P >1e-5) %>% group_by(SNP) %>% count() %>% filter(n==4) %>% select(SNP) %>% distinct() %>%
    write.table("hwe1e-5_superpops_195ksnps", sep="\t", row.names = F, quote = F)
    ' --args ${out}BED/
    #If we wanted to give the user options to change the MAF cutoff, we would need to expose this option in the filter
}

king_coefficients(){
    module load $plink2Load
    #HAVE I ALREADY MADE ${OUT}KING?
    mkdir -p ${out}KING/matrix
    plink2 \
    --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect \
    --extract ${out}BED/hwe1e-5_superpops_195ksnps \
    --make-king triangle bin \
    --out ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE1e_5 \
    --thread-num 30
    
    plink2 --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect \
    --king-cutoff ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE1e_5 0.0442 && \
    mv ${out}KING/plink2.king.cutoff.in.id  ${out}KING/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE1_5.king.cutoff.in.id && \
    mv ${out}KING/plink2.king.cutoff.out.id  ${out}KING/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE1_5.king.cutoff.out.id
    
    #Filter the file
    plink2 --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect \
    --make-bed \
    --keep ${out}KING/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE1e_5.king.cutoff.in.id \
    --out ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_unrelated
    
    #Also produce a related set
    #Filter the file
    plink2 --bfile ${out}BED/autosomes_LD_pruned_1kgp3Intersect \
    --make-bed \
    --remove ${out}KING/autosomes_LD_pruned_1kgp3Intersect_triangle_HWE1e_5.king.cutoff.in.id \
    --out ${out}KING/matrix/autosomes_LD_pruned_1kgp3Intersect_related
}

##
# Run Thanos PCs script - look for this
##
1kgp3_ancestries(){
    module load bio/GCTA/1.93.1_beta
    cd /re_gecip/BRS/thanos/ethnicities_aggV2
    file=1KGP3_30K_unrel_autosomes
    folder=/re_gecip/BRS/thanos/ethnicities_aggV2/
    path=$folder$file
    out=1KGP3_30K_unrel_autosomes
    gcta64 --bfile $path --make-grm-bin --thread-num 30 --out $file
    gcta64 --grm $out --pca 20  --out $out --thread-num 30
    gcta64 --bfile $path --pc-loading $out --out $out --thread-num 30
}

##
# Run the Ancestry inference script
##
infer_ancestry(){
    #TODO - sort out the hardcoded filepaths
    module load $RLoad
    mkdir -p ${out}Ancestries
    export R_LIBS="$(dirname `which R`)/../lib64/R/library" #Set me to mitigate library mismatch on helix
    Rscript infer_ancestry.R  \
    ${resources}/${kgp3_sample_table} \
    ${resources}/${super_pop_codes} \
    ${resources}/${kgp3_unrel} \
    "/re_gecip/BRS/thanos/ethnicities_aggV2/aggV2_ethn_200Ksites/1KGP3_intersectGEL_200Kset_perpopHWE1e-6_unrel_maf0.05both1K100K.eigenvec" \
    "/re_gecip/BRS/thanos/ethnicities_aggV2/aggV2_ethn_200Ksites/1KGP3_intersectGEL_200Kset_perpopHWE1e-6_unrel_maf0.05both1K100K_GELprojection.proj.eigenvec" \
    ${out}Ancestries
}

## END OF ANCESTRY INFERENCE/KING WORKFLOW ##


p_hwe(){
    module load $bcftoolsLoad
    module load $plinkLoad
    mkdir -p ${out}HWE
    #Take the unrelated set, split by ancestry
    ####
    # What is needed for this?
    # Check the dat and unrel vars in the R script below. the dat is inferred ancestries per sample.
    # Unrels is given by Thanos
    # Once these are changed the script will run fine as is
    # You can run this interactively as computational overheads are v low
    ####
    R -e 'library(data.table)
        library(dplyr)
        args <- commandArgs(trailingOnly = T)
        dat <- fread(args[1]) %>% as_tibble()
        unrel <- fread(args[2]) %>% as_tibble()
        dat <- dat %>% filter(mafList.pcsgel.Sample %in% unrel$IID)
        for(col in c("AFR","EUR","SAS","EAS")){dat[dat[col]>0.8,c("mafList.pcsgel.Sample",col)] %>%
        write.table(paste0(args[3],col,"_unrelated_pop.txt"), quote = F, row.names=F, sep = "\t")};
    ' --args ${out}Ancestries/predicted_ancestries.tsv ${out}UNRELATEDLIST ${out}OUTLOC
    
    for pop in AFR EUR SAS EAS; do
        
        echo -e "Calculating HWE on ${pop} samples..."
        plink --bfile ${POINTS TO LIST OF HIGH Q SNPS IN PLINK FORMAT} \
        --hardy midp \
        --keep ${out}/ancestry/${pop}_unrelated_pop.keep \
        --double-id \
        --allow-extra-chr \
        --out ${out}/HWE/${i}_$pop
    done
}


endAggregateAnnotation2(){
    ###############
    ##  Purpose: Annotate and make pass/fail. print subset of cols
    ##  Input: All the outputs of siteQC metrics
    ##  Output: A single text file containing annotated variants
    ###############
    module load $RLoad
    export final='TRUE'
    autosomes=`grep -v chrX ${bedfile} | grep -v chrY | grep -v chrX | grep -v chrM`
    export infile=`echo $autosomes | cut -d' ' -f ${LSB_JOBINDEX}`
    #Only autosomes
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
    
    #Some of the below will need to be changed when we finalise filepaths etc for data
    mkdir -p ${out}Annotation_final
    Rscript annotatePerChunk.R \
    ${i} \
    ${out}startfile/start_file_${i} \
    ${out}missing2/missing1_${i} \
    ${out}missing2/missing2_${i} \
    ${out}medianCoverageAll/medianCoverageAll${i} \
    ${out}medianCovNonMiss/medianCovNon_${i} \
    ${out}medianGQ/medianGQ_${i} \
    ${out}AB_hetAll/hetAll_${i} \
    ${out}AB_hetPass/hetPass_${i} \
    ${out}MendelErrSites/MendErr_${i}.lmendel \
    ${out}Annotation_final \
    ${resources}/N_samples \
    ${out}AC_counts/${i}_AC \
    ignore #unused arg in this argument
    #Compress and tabix files
    
    bgzip -f ${out}Annotation_final/BCFtools_site_metrics_${i}.txt
    tabix -f -s1 -b2 -e2 ${out}Annotation_final/BCFtools_site_metrics_${i}.txt.gz
    
}

make_header(){
    #Need to correct paths for this
    [[ ! -f all_seen_flags.txt ]] && awk '{print $1}' all_summary_filter.txt | sort | uniq  > all_seen_flags.txt
    #Start with header from input file
    #Remember that we need to use the header with removed samples
    bcftools view -h all_chunks_merged_norm_chr11_52570777_53572719.bcf -S \
    ${resources}78389_final_platekeys_agg_v9.txt --force-samples > orig.hdr
    #Clean out INFO lines,  these wil be repopulated
    grep -v \#\#INFO orig.hdr > orig.hdr2 && mv orig.hdr2 orig.hdr
    
    #Now we need to insert our header lines
    tee headerbuild.R << EOF
    #This script runs interactively, doesn't work otherwise, not sure why
        library(data.table);library(dplyr);library(magrittr);library(stringr);
        #Read in the header options
        toAdd <- fread("all_seen_flags.txt", header = F) %>% as_tibble() %>%
        filter(V1 != 'PASS');
        starter <- '##FILTER=<ID=';
        mid <- ',Description="';
        end <- '.">'
        #Construct descriptions. We won't be all that descriptive for the filter field
        #more so for the info field. Just state failures for filter
        #Probably should add what the failing values are. Maybe only include these for the single instance onces
        #to save on space
        singlecases <- setNames(
            as.list(
                c('Fails depth, median depth < 10',
                'Fails GQ, median GQ <15',
                'Fails missingness, missingness > 0.05',
                'Fails completeGTRatio, compelete sites <0.5',
                'Fails AB ratio, AB ratio of (het sites passing binomial distribution with p-value < 0.01 / all het sites) < 0.25',
                'Fails phwe_eur, site is out of HWE (<10e-6) for inferred unrelated inferred eur superpopulation'
                )),
                c('depth',
                'GQ',
                'missingness',
                'completeGTRatio',
                'ABratio',
                'phwe_eur')
            )

        singlecases <- paste0(starter, names(singlecases), mid, singlecases, end)


    #Now do all the other filter combinations
        toAdd %<>% mutate(V2 = case_when(str_count(V1, ':') == 1 ~ str_replace(V1, ':',' and '),
                               str_count(V1,':') > 1 ~ str_replace_all(V1,':',', '),
                                TRUE ~ V1))
        toAdd %<>% mutate(V2 = ifelse(str_count(V2,',') > 1, sub(".([^,]*)$", " and\\1", V2), V2))


        toAdd %<>% mutate(V3 = paste0(starter,V1, mid, 'Fails ', V2, end))

        #Now sort out the info
        infostart <- '##INFO=<ID='
        infomid <- ',Number=.,Type=Float,Description="'
        infoend <- '">'

        infos <- setNames(
            as.list(
                    c('Median depth (taken from the DP FORMAT field) of all samples. Used for filter flag depth.',
                      'Median depth (taken from the DP FORMAT field) from samples with non-missing genotypes.',
                      'Median genotype quality(taken from the GQ FORMAT field) from samples with non-missing genotypes. Used for filter flag GQ.',
                      "Ratio of fully missing genotypes (GT = './.' and FORMAT/DP = 0)",
                      'The ratio of complete GTs/total number of samples',
                      'For each het call, a binomial test is conducted for reads supporting the ref and alt alleles. AB ratio is the hets showing imbalance (p<0.01) divided by the total number of hets.',
                      'The number of Mendel Errors at this site, calculated from confirmed trios',
                      'HWE mid p-value in inferred unrelated inferred afr superpop.',
                      'HWE mid p-value in inferred unrelated inferred amr superpop.',
                      'HWE mid p-value in inferred unrelated inferred eas superpop.',
                      'HWE mid p-value in inferred unrelated inferred eur superpop.',
                      'HWE mid p-value in inferred unrelated inferred sas superpop.'
                      )),
                    c("medianDepthAll",
                    "medianDepthNonMiss",
                    "medianGQ",
                    "missingness",
                    "completeGTRatio",
                    "ABratio",
                    "MendelSite",
                    "phwe_afr",
                    "phwe_amr",
                    "phwe_eas",
                    "phwe_eur",
                    "phwe_sas"
                    )
                )
            #Now build this info full string
            infosout <- paste0(infostart, names(infos), infomid, infos, infoend )

    #Now print this out to file, and add it to the rest of the header
    names(singlecases <- NULL)
    d <- c(singlecases, toAdd$V3, infosout) %>% unlist() %>% as.data.frame()
    fwrite(d, 'additional_header.txt', quote = F, col.names = F)

EOF
    Rscript headerbuild.R
}

annotate_bcfs(){
    #Write this out to different directory now.
    #Version to generalise
    #First prepare annotation file
    module load $bcftoolsLoad
    
    autosomes=`grep -v chrX ${bedfile} | grep -v chrY | grep -v chrX | grep -v chrM`
    export infile=`echo $autosomes | cut -d' ' -f ${LSB_JOBINDEX}`
    #Only autosomes
    export i=`echo $infile | awk -F "_" 'sub(/.vcf.gz/,"",$NF) {print $(NF-2)"_"$(NF-1)"_"$NF}'`
    
    #Now use final sample list to annotate bcf
    bcftools view \
    -S /gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/additional_data/sample_list/aggV2_sampleIds_mpv10_78195.tsv \
    ${input}${infile} -Ou \
    --force-samples |\
    bcftools annotate  \
    -x FILTER,^INFO/OLD_MULTIALLELIC,^INFO/OLD_CLUMPED \
    -a  ${out}Annotation_final_corrected/BCFtools_site_metrics_${i}.txt.gz \
    -h additional_header.txt \
    --threads 30 \
    -c CHROM,POS,REF,ALT,missingness,medianDepthAll,medianDepthNonMiss,medianGQ,completeGTRatio,MendelSite,ABratio,phwe_afr,phwe_eur,phwe_eas,phwe_sas,FILTER,- | \
    bcftools +fill-tags \
    -o /gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data/gel_mainProgramme_aggV2_${i}.vcf.gz \
    -Oz \
    --threads 30 \
    -- -d -t AC,AC_Hom,AC_Het,AC_Hemi,AN
    
    bcftools index /gel_data_resources/main_programme/aggregation/aggregate_gVCF_strelka/aggV2/genomic_data/gel_mainProgramme_aggV2_${i}.vcf.gz
}


#Export functions

#Step 1 metric calculations
export -f startFile
export -f missingness1
export -f missingness2
export -f medianCovAll
export -f ABRatioP1
export -f ABRatioP2
export -f completeSites
export -f medianGQ
export -f medianCovNonMiss
export -f aggregateAnnotation
export -f pull_1KGPSites

#King workflow



#Step 2 metric calculations

#Aggregation

