import os, glob, sys, re

from snakemake.shell import shell
import logging, traceback
from pathlib import Path
from collections import defaultdict

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


# Install exception handler
sys.excepthook = handle_exception

# Start script
sample_table = snakemake.params.sample_table
parse_sample_table = snakemake.params.parse_sample_table
new_parse_sample_table = snakemake.output.new_parse_sample_table
fastq_info = snakemake.params.fastq_info
failed_sample = snakemake.params.failed_sample

def is_paired_fastq(file_path):
    pattern = re.compile(r'(.+?)_paired_([12])\.fastq(?:\.gz)?$')
    return bool(pattern.search(str(file_path)))

def get_parse_sample_table_by_sample_table(sample_table, failed_sample):
    if not Path(sample_table).exists(): 
        logging.error(f"{sample_table} does not exist!")
        sys.exit(1)
    all_sample2infos = {}
    with open(sample_table, "r") as f:
        for line in f:
            if line.startswith("#"): continue
            tokens = line.strip().split("\t")
            assert len(tokens) >= 3
            projects = tokens[1]
            if "," in projects:
                projects = projects.split(",")
            elif ";" in projects:
                projects = projects.split(";")
            else:
                projects = projects.split()
            
            reads_path = Path(tokens[2])
            reads_path_name = reads_path.name

            projects2num = {}
            if projects[0] in reads_path_name:
                reads_path_name_tokens = reads_path_name.split("_")
                assert len(reads_path_name_tokens) % 2 == 0
                for i in range(0, len(reads_path_name_tokens), 2):
                    projects2num[reads_path_name_tokens[i]] = int(reads_path_name_tokens[1])
            sample_num = sum(projects2num.values())
            sample2path = defaultdict(list)
            for file_path in reads_path.rglob("*"):
                if not is_paired_fastq(file_path): continue
                if file_path.is_file() and (str(file_path).endswith(("fastq.gz", "fastq", "fq.gz", "fq"))):
                    sample = file_path.stem.split("_")[0].split(".")[0]
                    # if sample in sample2path:
                    #     if sample2path[sample].endswith("gz"):
                    #         sample2path[sample] = file_path
                    #     else: 
                    #         continue
                    # else:
                    #     sample2path[sample] = file_path
                    sample2path[sample].append(str(file_path))
            # if sample_num != len(sample2path):
            #     logging.info(f"project {tokens[1]} sample number {len(sample2path)} is not equal to real number {sample_num}")
            #     print(f"project {tokens[1]} sample number {len(sample2path)} is not equal to real number {sample_num}")


            for sample, paths in sample2path.items():
                if Path(paths[0]).stem.endswith(("_1", "_2")): 
                    paired = True
                else: 
                    paired = False
                has_gz = any('.gz' in f for f in paths)
                if has_gz:
                    fq_files = [f for f in paths if f.endswith('.gz')]
                    if paired:
                        try:
                            assert len(fq_files) == 2
                        except:
                            print(fq_files)
                else:
                    fq_files = [f for f in paths if not f.endswith('.gz')]
                    if paired:
                        try:
                            assert len(fq_files) == 2
                        except:
                            print(fq_files)
                fq_files.sort()
                project_name = "_".join(projects)
                all_sample2infos[sample] = [project_name, sample, ",".join(fq_files), tokens[0], tokens[1]]
    with open(failed_sample, "r") as f:
        failed_sample_list = f.read().splitlines()
    if len(failed_sample_list) > 0:
        for sample in failed_sample_list:
            try:
                del all_sample2infos[sample]
            except KeyError as e:
                print(f"KeyError: {e}") 
    return all_sample2infos

def get_parse_sample_table_in_paras(parse_sample_table):
    all_sample2infos = {}
    if not Path(parse_sample_table).exists():
        return all_sample2infos
    with open(parse_sample_table, "r") as f:
        for line in f:
            if line.startswith("#"): continue
            tokens = line.strip().split("\t")
            all_sample2infos[tokens[1]] = tokens
    return all_sample2infos


new_all_sample2infos = get_parse_sample_table_by_sample_table(sample_table, failed_sample)
exist_all_sample2infos = get_parse_sample_table_in_paras(parse_sample_table)

if not Path(fastq_info).exists:
    Path(fastq_info).mkdir(parents=True)

if new_all_sample2infos != exist_all_sample2infos:
    need_update_samples = []
    
    for sample, new_infos in new_all_sample2infos.items():
        if sample not in exist_all_sample2infos:
            need_update_samples.append(sample)
        elif exist_all_sample2infos[sample] != new_infos:
            need_update_samples.append(sample)
    
    removed_samples = []
    for sample in exist_all_sample2infos:
        if sample not in new_all_sample2infos:
            removed_samples.append(sample)
    
    exist_all_sample2infos.update(new_all_sample2infos)
    
    for sample in removed_samples:
        if sample in exist_all_sample2infos:
            del exist_all_sample2infos[sample]
    
    with open(new_parse_sample_table, "w") as f:
        f.write("#project_name\taccession\tpath\tdisease\tproject\n")
        for sample, infos in sorted(exist_all_sample2infos.items()):
            f.write("\t".join(infos) + "\n")
    
    for sample in need_update_samples:
        if sample in exist_all_sample2infos:
            sample_path = Path(fastq_info) / sample
            sample_path.mkdir(parents=True, exist_ok=True)
            
            sample_file = sample_path / f"{sample}.txt"
            with open(sample_file, "w") as f2:
                infos = exist_all_sample2infos[sample]
                if len(infos) > 2:
                    file = infos[2].replace(",", "\n")
                    f2.write(file + "\n")



            



                                                
            

            

