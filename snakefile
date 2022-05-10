    
# alias smk='mv logs/*.txt logs/old 2> /dev/null; snakemake --profile profiles/slurm'
# ulog; snakemake --profile profiles/local/ --until msfragger -p                                                                    

__author__: "Carl Mathias Kobel & Arturo Vera De Ponce Leon"



# who wrote this


import os
from datetime import datetime
import time
import re
#from shutil import copyfile
import glob

import pandas as pd
import re


print("/*                                                                               ")
print("                            _________________________________________            ")
print("                         <               MS_pipeline2                 >          ")
print("                         < Taking what worked and leaving what didn't >          ")
print("                           ------------------------------------------            ")
print("                                                          \\                     ")
print("                             ___......__             _     \\                    ")
print("                         _.-'           ~-_       _.=a~~-_                       ")
print(" --=====-.-.-_----------~   .--.       _   -.__.-~ ( ___===>                     ")
print("               '''--...__  (    \\ \\\\\\ { )       _.-~                         ")
print("                         =_ ~_  \\\\-~~~//~~~~-=-~                               ")
print("                          |-=-~_ \\\\   \\\\                                     ")
print("                          |_/   =. )   ~}                                        ")
print("                          |}      ||                                             ")
print("                         //       ||                                             ")
print("                       _//        {{                                             ")
print("                    '='~'          \\\\_    =                                    ")
print("                                    ~~'                                          ")
print("                                                                                 ")


# TODO: It looks like there is a problem with annotate. Should it really be running for each sample individually? I would think that once per batch should be plentiful.



# Import configuration
configfile: "config.yaml"
config_batch = config["batch"]
config_d_base = config["batch_parameters"][config_batch]["d_base"]
config_database_glob = config["batch_parameters"][config_batch]["database_glob"]
config_database_glob_read = glob.glob(config_database_glob)
config_samples = config["batch_parameters"][config_batch]["samples"]

# Present configuration
print(f"config_batch:         '{config_batch}'")
print(f"config_d_base:        '{config_d_base}'")
print(f"config_database_glob: '{config_database_glob}:'")
for i, j in enumerate(config_database_glob_read):
    print(f"  {i+1}) {j}")
print()


# Populate dataframe
df = pd.DataFrame(data = {'sample':  config_samples.keys(),
                          'barcode': config_samples.values()})

df["basename"] = [re.sub(".d$", "", barcode) for barcode in df["barcode"]]
#df["path"] = config_d_base + "/" + df["barcode"]






print(df)
print("//")
print()




# Define default workflow
rule all:
    input: expand(["output/{config_batch}/metadata.tsv", \
                   "output/{config_batch}/{basename}.pepXML"], \
                   config_batch = config_batch, \
                   sample = df["sample"], \
                   basename = df["basename"])

#"output/{config_batch}/samples/{sample}/peptideprophet-{sample}.pep.xml", \
#                   "output/{config_batch}/samples/{sample}/protein.tsv", \
#                   "output/{config_batch}/samples/{sample}/{sample}_quant.csv"], \
#
# #"output/{config_batch}/samples/{sample}/annotate.done", \ # removed this one in the name of using the .meta/ files instead.

#"output/{config_batch}/samples/{sample}/{sample}_quant.csv"], \ # disabled until msfragger basename-sample conversion works

#"output/{config_batch}/msfragger/philosopher_database.fas", 






# Simply store the batch metadata in the output directory for later reference.
rule metadata:
    input: "output/{config_batch}/msfragger/link_input.done"
    output: "output/{config_batch}/metadata.tsv"
    params: dataframe = df.to_csv(None, index_label = "index", sep = "\t")
    shell: """

        echo '''{params.dataframe}''' > {output}
    
    """




# Instead of moving the output files from the backup directory where the raw data is stored, I'd rather link the raw data to the end destination and not have to move anything. Hence this rule.
rule link_input:
    input:
        d_files = directory((config_d_base + "/" + df["barcode"]).tolist())
    output:
        linked_d_files = directory("output/{config_batch}/" + df["barcode"]), # Not needed when we have the linked_flag file.
        linked_flag = touch("output/{config_batch}/link_input.done"),
    params:
        dir = "output/{config_batch}" # Why not make this as an actual output instead. # Because then you can't have outputs in this dir before this rule
    shell:"""

        ln -s {input.d_files} {params.dir}

        """




# Build a database of the known amino acid sequences from the selected database glob.
rule philosopher_database:
    input: glob.glob(config_database_glob)
    output: 
        database = "output/{config_batch}/philosopher_database.fas",
        #meta = "output/{config_batch}/msfragger/.meta/db.bin"
    #benchmark: "output/{config_batch}/benchmarks/database.tab"
    threads: 8
    params:
        philosopher = config["philosopher_executable"]
    shell: """


        >&2 echo "Catting database files ..."
        # Cat all database source files into one.
        cat {input} > output/{config_batch}/cat_database_sources.faa


        >&2 echo "Change dir ..."
        # As philosopher can't specify output files, we need to change dir.
        mkdir -p output/{config_batch}
        cd output/{config_batch}

        >&2 echo "Philosopher workspace clean ..."
        {params.philosopher} workspace \
            --nocheck \
            --clean 

        >&2 echo "Philosopher workspace init ..."
        {params.philosopher} workspace \
            --nocheck \
            --init 

        >&2 echo "Removing previous .fas ..."
        rm *.fas || echo "nothing to delete" # Remove all previous databases if any.

        >&2 echo "Philosopher database ..."
        {params.philosopher} database \
            --custom cat_database_sources.faa \
            --contam 




        >&2 echo "Move output ..."
        # Manually rename the philosopher output so we can grab it later
        mv *-decoys-contam-cat_database_sources.faa.fas philosopher_database.fas

        >&2 echo "Clean up ..."
        # clean up 
        rm cat_database_sources.faa


        """






# TODO: Test the implications of using shadow: "minimal"
rule msfragger:
    input:
        linked_flag = "output/{config_batch}/link_input.done",
        database = "output/{config_batch}/philosopher_database.fas",  
        d_files = ("output/{config_batch}/" + df["barcode"]).tolist()
    output:
        pepXMLs = "output/{config_batch}/" + df["basename"] + ".pepXML",
        stdout = "output/{config_batch}/msfragger.out.txt"
    # Using shadow is a nice way of getting rid of irrelevant outputs.
    #shadow: "minimal" # The setting shadow: "minimal" only symlinks the inputs to the rule. Once the rule successfully executes, the output file will be moved if necessary to the real path as indicated by output.
    # Shadow doesn't work well with tee, as tee needs access to the log directory
    threads: 8
    params:
        config_d_base = config_d_base,
        msfragger_jar = config["msfragger_jar"],
        n_samples = len(df.index)
    conda: "envs/openjdk.yaml"
    shell: """

        >&2 echo "MSFragger ..."
        java \
            -Xmx64G \
            -jar {params.msfragger_jar} \
            --num_threads {threads} \
            --database_name {input.database} \
            --output_location "output/{wildcards.config_batch}/" \
            {input.d_files} \
            | tee output/{wildcards.config_batch}/msfragger.out.txt

   

        # The last line extracts the number of lines that corresponds to the number of samples, to extract the total number of scans.
        # Without tee, the stdout would be lost.



        ls output/{wildcards.config_batch} > output/{wildcards.config_batch}/msfragger.done


        # makes a .pepindex and a pepXML for each sample.
        # I feel like it also creates a .mgf and .mzBIN in the source directory where the .d-dirs reside
        # Should I not move the .pepXML files? No, because I'm using the --output_location argument.
        # The tutorial mentions something about moving some .tsv files after running msfragger, but I haven't seen any.
        # These output files should in theory be mitigated by using the shadow rule?



        """


# Annotate for all files in one go
rule database_annotate:
    input: "output/{config_batch}/philosopher_database.fas",
    output: touch("output/{config_batch}/database_annotate.done")
    shell: """
        
        >&2 echo "Creating and changing dir ..."
        mkdir -p output/{config_batch}
        cd output/{config_batch}

        >&2 echo "Philosopher annotate ..."
        # Initially I tried preparing for all samples individually, but I will not do that anymore.
        {params.philosopher} database \
            --annotate ../../{input} 

    """


# I would like to run this per batch, but I'm afraid that the stupid output file naming will make too many conflics, so I don't dare to.
rule prophet_filter:
    input:
        database = "output/{config_batch}/msfragger/philosopher_database.fas",
        flag = "output/{config_batch}/database_annotate.done"
        #pepXML = lambda wildcards: "output/" + config_batch + "/msfragger/" + df[df["sample"] == wildcards.sample]["basename"] + ".pepXML",
        pepXMLs = "output/{config_batch}/" + df["basename"] + ".pepXML"
    #output: ["output/{config_batch}/samples/{sample}/ion.tsv", \
    #    "output/{config_batch}/samples/{sample}/peptideprophet-{sample}.pep.xml", \
    #    "output/{config_batch}/samples/{sample}/peptide.tsv", \
    #    "output/{config_batch}/samples/{sample}/protein.fas", \
    #    "output/{config_batch}/samples/{sample}/proteinprophet-{sample}.prot.xml", \
    #    "output/{config_batch}/samples/{sample}/protein.tsv", \
    #    "output/{config_batch}/samples/{sample}/psm.tsv"]
    output: flag("output/{config_batch}/prophet_filter_whatever.done")
    params:
        philosopher = config["philosopher_executable"]
    shell: """

    
        cd output/{wildcards.config_batch}/samples/{wildcards.sample}



        >&2 echo "Peptideprophet ..."
        {params.philosopher} peptideprophet \
            --nonparam \
            --expectscore \
            --decoyprobs \
            --ppm \
            --accmass \
            --output peptideprophet \
            --database ../../../../{input.database} \
            {wildcards.sample}.pepXML

        # Be noted that the warning about a missing file is not critical: https://github.com/cmkobel/MS_pipeline1/issues/2


FATAL: Because proteinprophet puts the output WHERE THE WORKSPACE WAS INITIALLY CREATED (prepare_annotate/), there will simply be too many hoops and lumps for this project to ever give fruit. Goodbye
        >&2 echo "Proteinprophet ..."
        {params.philosopher} proteinprophet \
            --output proteinprophet-{wildcards.sample} \
            peptideprophet-{wildcards.sample}.pep.xml


        >&2 echo "Filter ..." 
        {params.philosopher} filter \
            --sequential \
            --razor \
            --mapmods \
            --pepxml peptideprophet-{wildcards.sample}.pep.xml \
            --protxml proteinprophet-{wildcards.sample}.prot.xml

        # Assuming that philosopher filter works in place
        # TODO: Ask Arturo if that is true.
        >&2 echo "Report ..."
        {params.philosopher} report

        
    """



# TODO: This rule ought to output the abundances named in the samples name and not the basename? I don't really see any neat way to do that 
rule ionquant:
    input:
        irrelevant = ["output/{config_batch}/samples/{sample}/ion.tsv", \
            "output/{config_batch}/samples/{sample}/peptide.tsv", \
            "output/{config_batch}/samples/{sample}/protein.fas", \
            "output/{config_batch}/samples/{sample}/proteinprophet-{sample}.prot.xml", \
            "output/{config_batch}/samples/{sample}/protein.tsv"], 
        psm = "output/{config_batch}/samples/{sample}/psm.tsv",
        pepXML = lambda wildcards: "output/" + config_batch + "/msfragger/" + df[df["sample"] == wildcards.sample]["basename"] + ".pepXML",

    output: #touch("output/{config_batch}/samples/{sample}/ionquant.done")
        csv = "output/{config_batch}/samples/{sample}/{sample}_quant.csv"
    threads: 8
    conda: "envs/openjdk.yaml"
    params:
        ionquant_jar = config["ionquant_jar"],
        config_d_base = config_d_base, # I think this one is global, thus does not need to be params-linked.
        basename = lambda wildcards: df[df["sample"] == wildcards.sample]["basename"].values[0]



    shell: """


        >&2 echo "Ionquant ..."
        java \
            -Xmx32G \
            -jar {params.ionquant_jar} \
            --threads {threads} \
            --psm {input.psm} \
            --specdir {params.config_d_base} \
            {input.pepXML} 
            # address to msfragger pepXML file


            # TODO: Ask Arturo if it makes any sense that I'm not using the pepXML from peptideprophet, but the one directly from msfragger

            # Apparently, --specdir should point to the msfragger pepxmls. Maybe, I just need to point to the msfragger dir.
            # Or maybe I need to point directly to the file.
            # --specdir output/220315_test/msfragger/20220302_A1_Slot1-01_1_1592.pepXML 
            # Maybe the other pepxml is the culprit

            #mv output/{config_batch}/msfragger/{wildcards.sample}_quant.csv output/{config_batch}/samples/{wildcards.sample}/{wildcards.sample}_quant.csv
            mv output/{config_batch}/msfragger/{params.basename}_quant.csv output/{config_batch}/samples/{wildcards.sample}/{wildcards.sample}_quant.csv


    """

rule rmarkdown:
    input:
        metadata = "output/{config_batch}/metadata.tsv",
        psms = "output/{config_batch}/samples/{sample}/psm.tsv",
        quants = "output/{config_batch}/samples/{sample}/{sample}_quant.csv", # This simply makes it only run if rule ionquant was successful.
    output:
        "output/{config_batch}/QC.html"
    conda: "envs/r-markdown.yaml"
    shell: """
        #Rscript --what scripts/QC.Rmd {input.metadata}

        #cp scripts/QC.Rmd QC.Rmd

        Rscript -e 'library(rmarkdown); rmarkdown::render("scripts/QC.rmd", "html_document")'

        #rm rmarkdown_template.rmd
        #mv rmarkdown_template.html ../{output}

"""

        


print("*/") # This is a language specific comment close tag that helps when you export the workflow as a graph



# TODO: Go through the whole pipeline one job at a time, and make sure that all outputs are managed in the rules.
# TODO: Export the dag and put it into the readme with a bit of documentation.