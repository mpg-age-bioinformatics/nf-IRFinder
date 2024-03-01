#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """
    if [[ "${params.containers}" == "singularity" ]] ;

      then

        cd ${params.image_folder}

        if [[ ! -f irfinder-1.3.1.sif ]] ;
          then
            singularity pull irfinder-1.3.1.sif docker://index.docker.io/mpgagebioinformatics/irfinder:1.3.1  
        fi

        if [[ ! -f rnaseq.python-3.8-2.sif ]] ;
          then
            singularity pull rnaseq.python-3.8-2.sif docker://index.docker.io/mpgagebioinformatics/rnaseq.python:3.8-2
        fi

        if [[ ! -f deseq2-1.38.0.sif ]] ;
          then
            singularity pull deseq2-1.38.0.sif docker://index.docker.io/mpgagebioinformatics/deseq2:1.38.0
        fi

    fi

    if [[ "${params.containers}" == "docker" ]] ;

      then
        docker pull mpgagebioinformatics/irfinder:1.3.1
        docker pull mpgagebioinformatics/rnaseq.python:3.8-2
        docker pull mpgagebioinformatics/deseq2:1.38.0

    fi
    """

}

process copy_IRFinder {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.project_folder}/IRFinder-1.3.1/bin/IRFinder").exists() )

  script:
    """
    rsync -rtvhl /IRFinder-1.3.1 ${params.project_folder}
    """
}

process small_int_script {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val small_introns

  output:
    val small_introns

  when:
    ( "${small_introns}" == "YES")

  script:
    """
    #!/usr/local/bin/python
with open("/workdir/IRFinder-1.3.1/bin/util/IntronExclusion.pl", 'r+') as file:
  lines = file.readlines()
  file.seek(0)
  for line_num, line in enumerate(lines, start=1):
      
    if line_num == 82 and line == '  if (\$newlen > 40 && (\$newlen/\$len) >= 0.7) {\\n':
        modified_line = '  if (\$newlen > 10 && (\$newlen/\$len) >= 0.5) {\\n'
        file.write(modified_line)
    elif line_num == 94 and line == '    if (\$len >= 110) {\\n':
        modified_line = '    if (\$len >= 30) {\\n'
        file.write(modified_line)
    elif line_num == 95 and line == '      print OF50 join("\\\\t", \$chr, \$start+5, \$start+55, "S", 0, \$dir, \$start+5, \$start+55, "255,0,0", 1, 50, 0), "\\\\n";\\n':
        modified_line = '      print OF50 join("\\\\t", \$chr, \$start+5, \$start+15, "S", 0, \$dir, \$start+5, \$start+15, "255,0,0", 1, 10, 0), "\\\\n";\\n'
        file.write(modified_line)
    elif line_num == 96 and line == '      print OF50 join("\\\\t", \$chr, \$end-55, \$end-5, "S", 0, \$dir, \$end-55, \$end-5, "255,0,0", 1, 50, 0), "\\\\n";\\n':
        modified_line = '      print OF50 join("\\\\t", \$chr, \$end-15, \$end-5, "S", 0, \$dir, \$end-15, \$end-5, "255,0,0", 1, 10, 0), "\\\\n";\\n'
        file.write(modified_line)
    else:
        file.write(line)
    """
}


process build_ref {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """
    mkdir -p "${params.project_folder}/REF"
    cd "${params.project_folder}/REF"
    ln -s ${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.fa genome.fa
    ln -s ${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.gtf transcripts.gtf

    /workdir/IRFinder-1.3.1/bin/IRFinder -m BuildRefProcess -r "${params.project_folder}/REF"
    """
}


process quantify_ir {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    tuple val(pair_id), path(fastq)

  output:
    val pair_id

  when:
    ( ! file("${params.project_folder}/irquant_out/${pair_id}/IRFinder-IR-dir.txt").exists() ) 
  
  script:
  def single = fastq instanceof Path

  if ( single ) {
    """
    mkdir -p "${params.project_folder}irquant_out"
    #cd "${params.IRquant_raw_data}"
    cd "${params.project_folder}raw_data/"

    /workdir/IRFinder-1.3.1/bin/IRFinder -r "${params.project_folder}/REF" -d ${params.project_folder}irquant_out/${pair_id} ${pair_id}.READ_1.fastq.gz
    """
  }
  else {
    """
    mkdir -p "${params.project_folder}irquant_out"
    cd "${params.project_folder}raw_data/"

    /workdir/IRFinder-1.3.1/bin/IRFinder -r "${params.project_folder}/REF" -d ${params.project_folder}irquant_out/${pair_id} ${pair_id}.READ_1.fastq.gz ${pair_id}.READ_2.fastq.gz
    """
  }
}

process make_comps {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.scripts}group_comparison.txt").exists() ) 
  
  script:
    """
    #!/usr/local/bin/python
import pandas as pd
import os

samplestable = pd.read_excel("${params.scripts}sample_sheet.xlsx")

ir_folder="${params.project_folder}/irquant_out"
ir_diff="${params.project_folder}/irdiff_out"
sufix="${params.read1_sufix}"

# count replicates per group
replicates = samplestable["group"].value_counts()
# separate into, no replicates, 2 replicates, more than 3 replicates
REP1 = samplestable[samplestable["group"].isin(replicates.index[replicates == 1])]
REP2 = samplestable[samplestable["group"].isin(replicates.index[replicates == 2])]
REP3 = samplestable[samplestable["group"].isin(replicates.index[replicates >= 3])]

##### create comparison list for no replicates
##### run as for line in file:
##### bin/analysisWithNoReplicates.pl per line

O = open("${params.scripts}noreplicates.txt", "w")
if REP1.shape[0] >= 2:
    # write out all pairwise comparisons
    for i, sample_i in REP1.iterrows():
        for j, sample_j in REP1.iterrows():
            if i > j:
                samp1=sample_i["Files"].replace(sufix, "")
                samp2=sample_j["Files"].replace(sufix, "")
                g1=sample_i["group"]
                g2=sample_j["group"]
                if ( os.path.isfile(f"{ir_folder}/{samp1}/IRFinder-IR-dir.txt") ) & ( os.path.isfile(f"{ir_folder}/{samp2}/IRFinder-IR-dir.txt") ) :
                    O.write(f"-A {ir_folder}/{samp1}/IRFinder-IR-dir.txt -B {ir_folder}/{samp2}/IRFinder-IR-dir.txt > {ir_diff}/{g1}.vs.{g2}.txt\\n" )
                else:
                    O.write(f"-A {ir_folder}/{samp1}/IRFinder-IR-nondir.txt -B {ir_folder}/{samp2}/IRFinder-IR-nondir.txt > {ir_diff}/{g1}.vs.{g2}.txt\\n" )
O.close()

#### for now just handle them as if they had no replicates, in the future, change this
#### create comparison list for no replicates
#### run as for line in file:
#### bin/analysisWithNoReplicates.pl per line

O = open("${params.scripts}noreplicates.txt", "a")
if REP2.shape[0] >= 2:
    # write out all pairwise comparisons
    for i, sample_i in REP2.iterrows():
        for j, sample_j in REP2.iterrows():
            if i < j:
                samp1=sample_i["Files"].replace(sufix, "")
                samp2=sample_j["Files"].replace(sufix, "")
                g1=sample_i["Files"].replace(sufix, "")
                g2=sample_j["Files"].replace(sufix, "")
                if ( os.path.isfile(f"{ir_folder}/{samp1}/IRFinder-IR-dir.txt") ) & ( os.path.isfile(f"{ir_folder}/{samp2}/IRFinder-IR-dir.txt") ) :
                    O.write(f"-A {ir_folder}/{samp1}/IRFinder-IR-dir.txt -B {ir_folder}/{samp2}/IRFinder-IR-dir.txt > {ir_diff}/{g1}.vs.{g2}.txt\\n" )
                else:
                    O.write(f"-A {ir_folder}/{samp1}/IRFinder-IR-nondir.txt -B {ir_folder}/{samp2}/IRFinder-IR-nondir.txt > {ir_diff}/{g1}.vs.{g2}.txt\\n" )            
O.close()

# case of at least 3 replicates
# submit one run per line in group_comparisons.txt file

groups = replicates.index[replicates >= 3]
comparisons = open("${params.scripts}group_comparison.txt", "w")
for i, g1 in enumerate(groups):
    for j, g2 in enumerate(groups):
        if i < j:
            subsamples = samplestable[samplestable["group"].isin([g1, g2])]
            # write file list
            filelist = open(ir_diff+"/%s.vs.%s.filePaths.txt" %(g1, g2), "w")
            for k, row in subsamples.iterrows():
                samp=row["Files"].replace(sufix, "")
                if os.path.isfile(f"{ir_folder}/{samp}/IRFinder-IR-dir.txt") :
                    filelist.write(f"{ir_folder}/{samp}/IRFinder-IR-dir.txt\\n")
                else:
                    filelist.write(f"{ir_folder}/{samp}/IRFinder-IR-nondir.txt\\n")
            filelist.close()
            
            # write metafile
            subsamples["Files"] = subsamples["Files"].apply(lambda x: str(x).replace(sufix, ""))
            subsamples.rename({"Files":"SampleNames", "group": "Condition"}, axis = 'columns', inplace = True)
            subsamples.to_csv(ir_diff+"/%s.vs.%s.experiment.txt" %(g1, g2), sep = "\\t", index = False)
            
            # add comparison to comparison file
            comparisons.write("%s.vs.%s\\n" %(g1,g2))
            
comparisons.close()
    """
}

process noRep_diff {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val comp

  output:
    val comp

  when:
    comp.trim().length() > 0

  script:
  """
  /workdir/IRFinder-1.3.1/bin/analysisWithNoReplicates.pl ${comp}
  """
}

process wRep_diff {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val comp

  output:
    val comp

  when:
    ( ! file("${params.project_folder}tmp/${comp}.Rdata").exists() ) 

  script:
  """
  #!/usr/bin/Rscript
  print(Sys.getenv("R_LIBS_USER"))
  library(DESeq2)
  library(openxlsx)

  source("/workdir/IRFinder-1.3.1/bin/DESeq2Constructor.R")  #Load IRFinder-related function

  g1 = unlist(strsplit("${comp}", "[.]vs[.]"))[1]
  g2 = unlist(strsplit("${comp}", "[.]vs[.]"))[2]

  results = read.table("${params.project_folder}/irdiff_out/${comp}.filePaths.txt")
  paths = as.vector(results[, 'V1'])    # File names must be saved in a vector
  
  experiment = read.table("${params.project_folder}/irdiff_out/${comp}.experiment.txt",header=T)                       
  experiment[, 'Condition']=factor(experiment[, 'Condition'],levels=c(g1,g2))    # Set WT as the baseline in the analysis
  rownames(experiment)=NULL                         # Force removing rownames
  
  metaList=DESeqDataSetFromIRFinder(filePaths=paths, designMatrix=experiment, designFormula=~1)
  dds = metaList[['DESeq2Object']]                  # Extract DESeq2 Object with normalization factors ready
  
  design(dds) = ~Condition + Condition:IRFinder     # Build a formula of GLM. Read below for more details. 
  dds = DESeq(dds)                                  # Estimate parameters and fit to model
  
  # results for g1
  res.WT = as.data.frame(results(dds, name = paste0("Condition", g1, ".IRFinderIR")))
  res.WT[, "IR_vs_Splice"]=2^res.WT[,'log2FoldChange']
  res.WT[, 'IRratio'] = res.WT[, "IR_vs_Splice"]/(1+res.WT[, "IR_vs_Splice"])
  names(res.WT) = paste(g1, names(res.WT), sep ='.')
  
  #results for g2
  res.KO = as.data.frame(results(dds, name = paste0("Condition", g2,".IRFinderIR")))
  res.KO[, "IR_vs_Splice"] =2^res.KO[, 'log2FoldChange']
  res.KO[, "IRratio"] = res.KO[, "IR_vs_Splice"]/(1+res.KO[, "IR_vs_Splice"]) # That is weird, is that correct??
  names(res.KO) = paste(g2, names(res.KO), sep ='.')
  
  # Finally we can test the difference of (intronic reads/normal spliced reads) ratio between WT and KO
  res.diff = as.data.frame(results(dds, contrast=list(paste0("Condition", g2,".IRFinderIR"),paste0("Condition", g1, ".IRFinderIR"))))
  
  # We can plot the changes of IR ratio with p values
  # In this example we defined significant IR changes as
  # 1) IR changes no less than 10% (both direction) and 
  # 2) with adjusted p values less than 0.05
  
  res.diff[, 'IR.change'] = res.KO[, paste0(g2, ".IRratio")] - res.WT[, paste0(g1, '.IRratio')]
  names(res.diff) = paste0('diff.', names(res.diff))
  
  pdf("${params.project_folder}/irdiff_out/${comp}.pdf")
  plot(res.diff[, 'diff.IR.change'],col=ifelse(res.diff[, 'diff.padj'] < 0.05 & abs(res.diff[, 'diff.IR.change'])>=0.1, "red", "gray"), pch = 19, cex = 0.8, ylab = "IR.change")
  legend("bottomright", legend =  c('p < 0.05 and IRchange > 10%', 'n.s.'), pch = 19, col = c('red', 'gray'), bty = "n")
  
  plot(x = res.diff[, 'diff.log2FoldChange'], y = -log10(res.diff[, 'diff.padj']) ,col=ifelse(res.diff[, 'diff.padj'] < 0.05 & abs(res.diff[, 'diff.IR.change'])>=0.1, "red", "gray"), pch = 19, cex = 0.8, xlab = "IR.change", ylab = "-log10 padj")
  legend("top", legend =  c('p < 0.05 and IRchange > 10%', 'n.s.'), pch = 19, col = c('red', 'gray'), bty = "n")
  dev.off()
  
  # save tables
  res = cbind(res.diff, res.WT, res.KO)
  write.xlsx(res, "${params.project_folder}/irdiff_out/${comp}.xlsx", row.names = TRUE)
  
  save.image("${params.project_folder}/tmp/${comp}.Rdata")
  
  sink('${params.project_folder}/tmp/${comp}.sessionInfo.txt')
  sessionInfo()
  sink()
  """
  
}

process upload_paths {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
  """
    rm -rf upload.txt

    cd ${params.project_folder}/irdiff_out

    for f in \$(ls *) ; do echo "irdiff \$(readlink -f \${f})" >>  upload.txt_ ; done

    uniq upload.txt_ upload.txt 
    rm upload.txt_
    
    cd ${params.project_folder}/irquant_out

    for d in */ ; do  
      echo \${d%/}
      tar --exclude=\${d%/}/Unsorted.bam -cvzf \${d%/}.tar.gz \$d*
      echo "irquant_\${d%/} \$(readlink -f  \${d%/}.tar.gz)" >> upload.txt_
    done

    uniq upload.txt_ upload.txt 
    rm upload.txt_
  """
}

workflow images {
  main:
    get_images()
}

workflow upload {
  main:
    upload_paths()
}

workflow repo_IRFinder {
  main:
    copy_IRFinder()
}

workflow run_small_intron_calc {
  main:
    small_int_script( ${params.small_intron_calc} )
}

workflow run_build_ref {
  main:
    build_ref( )
}

workflow run_quantify_ir {
  main:
    read_files=Channel.fromFilePairs( "${params.project_folder}raw_data/*.READ_{1,2}.fastq.gz", size: -1 )
    quantify_ir( read_files )
    
}

workflow run_make_comps {

  if ( ! file("${params.project_folder}irdiff_out").isDirectory() ) {
    file("${params.project_folder}irdiff_out").mkdirs()
  }

  make_comps()
}

workflow run_noRep_diff {
  comps=Channel.fromPath("${params.scripts}/noreplicates.txt").splitText().map { it.trim() }
  noRep_diff( comps )
}

workflow run_wRep_diff {
  if ( ! file("${params.project_folder}tmp").isDirectory() ) {
    file("${params.project_folder}tmp").mkdirs()
  }

  comps=Channel.fromPath("${params.scripts}/group_comparison.txt").splitText().map { it.trim() }
  wRep_diff( comps )
}



