# MS_pipeline2



Fragpipe reminds me of a lego set I got when I was a kid, that did not come with a manual on how to assemble it. Here I'm trying to make my own manual in the form of a snakemake pipeline.


```
                            _________________________________________  
                         <               MS_pipeline2                 >
                         < Taking what worked and leaving what didn't >
                           ------------------------------------------  
                                                          \\           
                             ___......__             _     \\          
                         _.-'           ~-_       _.=a~~-_             
 --=====-.-.-_----------~   .--.       _   -.__.-~ ( ___===>           
               '''--...__  (    \\ \\\\\\ { )       _.-~               
                         =_ ~_  \\\\-~~~//~~~~-=-~                     
                          |-=-~_ \\\\   \\\\                           
                          |_/   =. )   ~}                              
                          |}      ||                                   
                         //       ||                                   
                       _//        {{                                   
                    '='~'          \\\\_    =                          
                                    ~~'                                
                                                                       
```




This is a very simple pipeline to transform .d files with an appropriate database into sane protein calls with abundances.

This pipeline can be seen as an enhanced version of https://fragpipe.nesvilab.org/docs/tutorial_linux.html

The workflow is based on msfragger and aims to keep a tidy tree of output files.


## Installation

1) Prerequisites:

  - A conda/miniconda3 environment
  - A conda environment with snakemake installed. You can install snakemake by following the [instructions in the official documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

2) Clone this repo

3) Set up the profiles/slurm/ configuration so that it matches your execution environment. There is also a profile for local execution without a job management system (profiles/local/)



## Usage

#### 1) Update config.yaml

The file config.yaml contains all the parameters needed for this pipeline to run successfully. You should update the parameters so they reflect which samples you wish to process.

Because nesvilab doesn't make their executables easily publicly available, you need to tell the pipeline where to find them on your system. Update addresses for the keys `philosopher_executable`, `msfragger_jar` and `ionquant_jar` which can be downloaded [here](https://github.com/nesvilab/philosopher/releases/latest), [here](https://github.com/Nesvilab/MSFragger/wiki/Preparing-MSFragger#Downloading-MSFragger) and [here](https://github.com/Nesvilab/IonQuant#download), respectively. 


Currently the pipeline only supports the input of .d-files ([agilent/bruker](https://en.wikipedia.org/wiki/Mass_spectrometry_data_format#Proprietary_formats)). Create an item in batch_parameters where you define key `d_base` which is the base directory where all .d-files reside. Define key `database_glob` which is a path (or glob) to the fasta-amino acid files that you want to include in the target protein database.

Define items under the `samples` key which link sample names to the .d-files.

Lastly, set the root `batch` key to point at the batch that you want to run.

#### 2) Run

Finally, run the pipeline in your command line with:
```
$ snakemake --profile profiles/slurm/ 
```

Below is an example graph for a run involving two samples.

![dag_](https://user-images.githubusercontent.com/5913696/167637412-48bebe40-a4a4-421c-a622-76e8f8470b3e.png)

## Future

This pipeline might involve an R-markdown performing trivial QC.





