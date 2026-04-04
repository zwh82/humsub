#!/usr/bin/env python3
"""
HuMSub Sourmash Analysis Pipeline
Distributed pipeline for subspecies classification using sourmash
"""

import os
import re
import sys
import yaml
import json
import argparse
import subprocess
import pandas as pd
from pathlib import Path
from datetime import datetime
import shlex
import concurrent.futures
import multiprocessing
from typing import List, Callable, Any, Tuple
import time
from collections import defaultdict


def load_config(config_path):
    """Load configuration from YAML file"""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config

def setup_directories(config):
    """Create necessary output directories"""
    output = config['output_dir']
    directories = [
        f"{output}/sketch_samples",
        f"{output}/gather",
        f"{output}/log/sketch_reads",
        f"{output}/log/gather",
    ]
    
    for directory in directories:
        Path(directory).mkdir(parents=True, exist_ok=True)
    
    return True

def download_sbt_index(config):
    """Download HuMSub SBT index"""
    sbt_path = os.path.join(config['resource_dir'], "HuMSub_51_1000.sbt.zip")
    
    if os.path.exists(sbt_path):
        print(f"SBT index already exists: {sbt_path}")
        return True
    
    print(f"Downloading SBT index to: {sbt_path}")
    cmd = f"wget -O {sbt_path} https://zenodo.org/records/15862096/files/HuMSub_51_1000.sbt.zip"
    
    try:
        subprocess.run(cmd, shell=True, check=True)
        print("Download completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Failed to download SBT index: {e}")
        return False

def get_samples_from_table(sample_table_path):
    """Extract sample names from sample table"""
    samples = []
    try:
        with open(sample_table_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                tokens = line.strip().split()
                if len(tokens) > 1:
                    samples.append(tokens[1])
    except Exception as e:
        print(f"Error reading sample table: {e}")
    
    return samples



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
    if Path(failed_sample).exists():
        with open(failed_sample, "r") as f:
            failed_sample_list = f.read().splitlines()
        if len(failed_sample_list) > 0:
            for sample in failed_sample_list:
                # try:
                #     del all_sample2infos[sample]
                # except KeyError as e:
                #     print(f"KeyError: {e}") 
                if sample in all_sample2infos:
                    del all_sample2infos[sample]
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

def check(config):
    sample_table = config["sample_table"]
    parse_sample_table = config["parse_sample_table"]
    fastq_info = config["fastq_info"]
    failed_sample = config["failed_sample"]
    new_parse_sample_table = config["parse_sample_table"]
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


def run_sketch_for_sample(sample, config):
    """Run sourmash sketch for a single sample"""
    fastq_dir = config["fastq_info"]
    output_dir = config['output_dir']
    reads_list = f"{fastq_dir}/{sample}/{sample}.txt"
    sketch_output = f"{output_dir}/sketch_samples/{sample}.sig"
    log_file = f"{output_dir}/log/sketch_reads/{sample}.log"
    
    if Path(sketch_output).exists():
        return True

    if not os.path.exists(reads_list):
        print(f"Error: Reads list not found for sample {sample}: {reads_list}")
        return False
    
    cmd = (
        f"sourmash sketch dna "
        f"-p k={config['kmer_len']},abund,scaled={config['scaled']} "
        f"--from-file {reads_list} "
        f"--merge {sample} "
        f"-o {sketch_output}"
    )
    
    try:
        with open(log_file, 'w') as log_f:
            result = subprocess.run(
                cmd,
                shell=True,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                text=True,
                check=False
            )
        
        success = result.returncode == 0
        print(f"Sample {sample}: {'success' if success else 'failed'}")
        return success
    except Exception as e:
        print(f"Sample {sample}: Exception - {e}")
        return False

def run_gather_for_sample(sample, config):
    """Run sourmash gather for a single sample"""
    output = config['output_dir']
    sketch_input = f"{output}/sketch_samples/{sample}.sig"
    gather_output = f"{output}/gather/{sample}.csv"
    log_file = f"{output}/log/gather/{sample}.log"
    
    if not os.path.exists(sketch_input):
        print(f"Error: Sketch file not found for sample {sample}: {sketch_input}")
        return False
    
    if Path(gather_output).exists():
        return True
    
    cmd = (
        f"sourmash gather "
        f"-k {config['kmer_len']} "
        f"--threshold-bp={config['threshold_bp']} "
        f"--scaled {config['scaled_downsample']} "
        f"-o {gather_output} "
        f"{sketch_input} {config['humsub_path']}"
    )
        
    try:
        with open(log_file, 'w') as log_f:
            result = subprocess.run(
                cmd,
                shell=True,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                text=True,
                check=False
            )
        
        success = result.returncode == 0
        print(f"Sample {sample}: {'success' if success else 'failed'}")
        return success
    except Exception as e:
        print(f"Sample {sample}: Exception - {e}")
        return False

def run_sketch_and_gather(sample, config):
    """Run sourmash sketch + gather for a single sample"""

    fastq_dir = config["fastq_info"]
    output_dir = config["output_dir"]

    # -------- paths --------
    reads_list = f"{fastq_dir}/{sample}/{sample}.txt"

    sketch_output = f"{output_dir}/sketch_samples/{sample}.sig"
    sketch_log = f"{output_dir}/log/sketch_reads/{sample}.log"

    gather_output = f"{output_dir}/gather/{sample}.csv"
    gather_log = f"{output_dir}/log/gather/{sample}.log"

    # -------- step 1: sketch --------
    if not Path(sketch_output).exists():
        if not os.path.exists(reads_list):
            print(f"Error: Reads list not found for sample {sample}: {reads_list}")
            return False

        sketch_cmd = (
            f"sourmash sketch dna "
            f"-p k={config['kmer_len']},abund,scaled={config['scaled']} "
            f"--from-file {reads_list} "
            f"--merge {sample} "
            f"-o {sketch_output}"
        )

        try:
            with open(sketch_log, 'w') as log_f:
                result = subprocess.run(
                    sketch_cmd,
                    shell=True,
                    stdout=log_f,
                    stderr=subprocess.STDOUT,
                    text=True
                )
            if result.returncode != 0:
                print(f"Sample {sample}: sketch failed")
                return False
            print(f"Sample {sample}: sketch success")
        except Exception as e:
            print(f"Sample {sample}: sketch exception - {e}")
            return False
    else:
        print(f"Sample {sample}: sketch skipped")

    # -------- step 2: gather --------
    if not Path(gather_output).exists():
        if not os.path.exists(sketch_output):
            print(f"Error: Sketch file missing for sample {sample}")
            return False

        gather_cmd = (
            f"sourmash gather "
            f"-k {config['kmer_len']} "
            f"--threshold-bp={config['threshold_bp']} "
            f"--scaled {config['scaled_downsample']} "
            f"-o {gather_output} "
            f"{sketch_output} {config['humsub_path']}"
        )

        try:
            with open(gather_log, 'w') as log_f:
                result = subprocess.run(
                    gather_cmd,
                    shell=True,
                    stdout=log_f,
                    stderr=subprocess.STDOUT,
                    text=True
                )
            if result.returncode != 0:
                print(f"Sample {sample}: gather failed")
                return False
            print(f"Sample {sample}: gather success")
        except Exception as e:
            print(f"Sample {sample}: gather exception - {e}")
            return False
    else:
        print(f"Sample {sample}: gather skipped")

    try:
        os.remove(sketch_output)
        print(f"Sample {sample}: sketch file removed")
    except Exception as e:
        print(f"Sample {sample}: failed to remove sketch file - {e}")

    return True

def run_gather_all(config):
    """Combine all gather results into a single file"""
    print("Combining all gather results...")
    
    # script_path = os.path.join(os.path.dirname(__file__), "scripts", "gather_all.py")
    script_path = Path(__file__).resolve().parent / "gather_all.py"
    if not os.path.exists(script_path):
        print(f"Error: gather_all.py script not found: {script_path}")
        return False
    
    cmd = f"python {str(script_path)} {config['output_dir']} {config['output_dir']}/subspecies_relab.csv {config['check_table']}"
    log_file = f"{config['output_dir']}/log/gather/all_gather.log"
    try:
        with open(log_file, 'w') as log_f:
            result = subprocess.run(
                cmd,
                shell=True,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                text=True,
                check=False
            )
        
        if result.returncode == 0:
            print("Successfully combined gather results")
            return True
        else:
            print(f"Failed to combine gather results: {result.stderr}")
            return False
    except Exception as e:
        print(f"Error running gather_all.py: {e}")
        return False

def run_map_taxonomy(config):
    """Map taxonomy to subspecies classification"""
    print("Mapping taxonomy...")
    
    # script_path = os.path.join(os.path.dirname(__file__), "scripts", "map_taxonomy.py")
    script_path = Path(__file__).resolve().parent / "map_taxonomy.py"
    
    if not os.path.exists(script_path):
        print(f"Error: map_taxonomy.py script not found: {script_path}")
        return False
    
    cmd = f"python {script_path} {config['output_dir']}/subspecies_relab.csv {config['taxonomy_table']} {config['taxonomy_version']} {config['output_dir']}/subspecies_taxonomy.csv"
    log_file = f"{config['output_dir']}/log/gather/map_taxonomy.log"
    try:
        with open(log_file, 'w') as log_f:
            result = subprocess.run(
                cmd,
                shell=True,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                text=True,
                check=False
            )
        
        if result.returncode == 0:
            print("Successfully mapped taxonomy")
            return True
        else:
            print(f"Failed to map taxonomy: {result.stderr}")
            return False
    except Exception as e:
        print(f"Error running map_taxonomy.py: {e}")
        return False

def process_samples_distributed(samples, config, task, threads, start_idx=None, end_idx=None):
    """
    Process a subset of samples for distributed execution
    
    Parameters:
    - samples: List of all samples
    - config: Configuration dictionary
    - task: 'sketch' or 'gather'
    - start_idx: Starting index (0-based)
    - end_idx: Ending index (exclusive)
    """
    if start_idx is None:
        start_idx = 0
    if end_idx is None:
        end_idx = len(samples)
    
    subset_samples = samples[start_idx:end_idx]
    
    print(f"Processing {len(subset_samples)} samples (indices {start_idx}-{end_idx-1})")
    success_count = 0
    if task == 'sketch':
        sketch_success = parallel_process(
            items=subset_samples,
            process_func=sketch_wrapper,
            func_name="sketch",
            max_workers=threads,
            config=config
        )
        print(f"Sketch: {sketch_success}/{len(subset_samples)} successful")
        success_count += sketch_success
    elif task == 'gather':
        gather_success = parallel_process(
            items=subset_samples,
            process_func=gather_wrapper,
            func_name="gather",
            max_workers=threads,
            config=config
        )
        print(f"Gather: {gather_success}/{len(subset_samples)} successful")
        success_count += gather_success
    elif task == "skga":
        gather_success = parallel_process(
            items=subset_samples,
            process_func=skga_wrapper,
            func_name="skga",
            max_workers=threads,
            config=config
        )
        print(f"Sketch and gather: {gather_success}/{len(subset_samples)} successful")
        success_count += gather_success        
    # success_count = 0
    # for i, sample in enumerate(subset_samples):
    #     print(f"[{i+1}/{len(subset_samples)}] Processing sample: {sample}")
        
    #     if task == 'sketch':
    #         success = run_sketch_for_sample(sample, config)
    #     elif task == 'gather':
    #         success = run_gather_for_sample(sample, config)
    #     else:
    #         print(f"Unknown task: {task}")
    #         success = False
        
    #     if success:
    #         success_count += 1
    
    # print(f"Completed {task} for {success_count}/{len(subset_samples)} samples")
    return success_count

def wrapped_func(args):
    process_func, item, func_kwargs = args
    try:
        return process_func(item, **func_kwargs)
    except Exception as e:
        return False, f"Unhandled exception: {str(e)}"

def parallel_process(
    items: List[Any],
    process_func: Callable,
    func_name: str = "process",
    max_workers: int = None,
    **func_kwargs
) -> Tuple[int, List[Tuple[Any, bool, str]]]:
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()
    
    total_items = len(items)
    print(f"Starting parallel {func_name} for {total_items} items")
    print(f"Using {max_workers} processes")
    
    start_time = time.time()
    success_count = 0
    
    # with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
    #     future_to_item = {
    #         executor.submit(wrapped_func, item): item
    #         for item in items
    #     }

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_item = {
            executor.submit(
                wrapped_func,
                (process_func, item, func_kwargs)
            ): item
            for item in items
        }

        processed_count = 0
        for future in concurrent.futures.as_completed(future_to_item):
            item = future_to_item[future]
            processed_count += 1
            
            percentage = (processed_count / total_items) * 100
            
            try:
                success = future.result()
                if success:
                    success_count += 1
            except Exception as e:
                message = f"Future exception: {str(e)}"
                print(message)
            
            elapsed = time.time() - start_time
            avg_time = elapsed / processed_count
            eta = avg_time * (total_items - processed_count)
            
            print(f"[{processed_count:3d}/{total_items:3d}] "
                  f"({percentage:5.1f}%) "
                  f"{str(item)[:30]:30s} "
                  f"| Elapsed: {elapsed:5.1f}s, ETA: {eta:5.1f}s")
    
    total_time = time.time() - start_time
    success_rate = (success_count / total_items) * 100
    print(f"\n{'='*80}")
    print(f"{func_name.upper()} COMPLETED")
    print(f"{'='*80}")
    print(f"Total items:        {total_items:>10}")
    print(f"Success:            {success_count:>10} ({success_rate:5.1f}%)")
    print(f"Failed:             {total_items - success_count:>10}")
    print(f"Total time:         {total_time:>10.1f} seconds")
    print(f"Average time/item:  {total_time/total_items:>10.2f} seconds")
    print(f"Processing rate:    {total_items/total_time:>10.2f} items/second")
    print(f"{'='*80}")
    
    return success_count

def sketch_wrapper(sample: str, config: dict) -> Tuple[bool, str]:
    print(f"sketch {sample}")
    return run_sketch_for_sample(sample, config)

def gather_wrapper(sample: str, config: dict) -> Tuple[bool, str]:
    return run_gather_for_sample(sample, config)

def skga_wrapper(sample: str, config: dict) -> Tuple[bool, str]:
    return run_sketch_and_gather(sample, config)

def main():
    print("Command:", " ".join(sys.argv))
    parser = argparse.ArgumentParser(description='HuMSub Sourmash Analysis Pipeline')
    parser.add_argument('--config', default='config/default_config.yaml', 
                       help='Path to configuration file')
    parser.add_argument('--task', choices=['all', 'download', 'process', 'sketch', 'gather', "skga", 'combine', 'taxonomy'],
                       default='all', help='Task to execute')
    parser.add_argument('--samples', nargs='+', help='Specific samples to process')
    parser.add_argument('--start', type=int, help='Start index for distributed processing')
    parser.add_argument('--end', type=int, help='End index for distributed processing')
    parser.add_argument('--skip-download', action='store_true', 
                       help='Skip SBT index download')
    parser.add_argument('--output-dir', default='output',
                       help='Output directory')
    parser.add_argument('--threads', default=64, type=int,
                       help='threads')     
    parser.add_argument('--check', default=None, type=str,
                       help='Path to configuration file')  
    parser.add_argument('--no-run', action='store_true',
                       help='Only check sample number')                         
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Set output directory in config
    config['output_dir'] = args.output_dir
    os.makedirs(config['output_dir'], exist_ok=True)
    output_dir = config['output_dir']
    # Set up directories
    setup_directories(config)
    
    # Set resource directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if not config.get('resource_dir', False):
        resource_dir = os.path.join(script_dir, "..", "resources")
        os.makedirs(resource_dir, exist_ok=True)
        config['resource_dir'] = resource_dir

    # Download SBT index if needed
    if not args.skip_download:
        if not download_sbt_index(config):
            print("Failed to download SBT index. Exiting.")
            sys.exit(1)
        if args.task == "download":
            sys.exit(0)
    
    check(config)
    # Parse samples
    samples = get_samples_from_table(config["parse_sample_table"])
    
    if not samples:
        print("No samples found. Exiting.")
        sys.exit(1)
    
    print(f"Found {len(samples)} samples")
    if args.no_run:
        sys.exit(0)
    # Filter samples if specific samples provided
    if args.samples:
        samples = [s for s in samples if s in args.samples]
        print(f"Processing {len(samples)} specific samples")
    

    # skga and process all perform sketch and gather, but the skga perform two steps on single sample at the same time.
    if args.task in ["all", "skga"]:
        print("\n" + "="*50)
        print("Running sketch and gather steps")
        Path(f"{output_dir}/sketch_samples").mkdir(parents=True, exist_ok=True)     
        Path(f"{output_dir}/gather").mkdir(parents=True, exist_ok=True)      
        if args.start is not None or args.end is not None:
            process_samples_distributed(samples, config, 'skga', args.threads, args.start, args.end)
            gather_success = parallel_process(
                items=samples,
                process_func=skga_wrapper,
                func_name="skga",
                max_workers=args.threads,
                config=config
            )
            print(f"Sketch and gather: {gather_success}/{len(samples)} successful")

    else:
        # Execute tasks
        if args.task in ['all', 'sketch']:
            print("\n" + "="*50)
            print("Running sketch step")
            print("="*50)
            # Create output directory
            Path(f"{output_dir}/sketch_samples").mkdir(parents=True, exist_ok=True)        
            if args.start is not None or args.end is not None:
                process_samples_distributed(samples, config, 'sketch', args.threads, args.start, args.end)
            else:
                # for sample in samples:
                #     run_sketch_for_sample(sample, config)
                sketch_success = parallel_process(
                    items=samples,
                    process_func=sketch_wrapper,
                    func_name="sketch",
                    max_workers=args.threads,
                    config=config
                )
                print(f"Sketch: {sketch_success}/{len(samples)} successful")

        if args.task in ['all', 'gather']:
            print("\n" + "="*50)
            print("Running gather step")
            print("="*50)
            # Create output directory
            Path(f"{output_dir}/gather").mkdir(parents=True, exist_ok=True)        
            if args.start is not None or args.end is not None:
                process_samples_distributed(samples, config, 'gather', args.start, args.end)
            else:
                # for sample in samples:
                #     run_gather_for_sample(sample, config)
                gather_success = parallel_process(
                    items=samples,
                    process_func=gather_wrapper,
                    func_name="gather",
                    max_workers=args.threads,
                    config=config
                )
                print(f"Gather: {gather_success}/{len(samples)} successful")

        if args.task in ['all', 'process']:
            print("\n" + "="*50)
            print("Running sketch and gather steps")
            print("="*50)
            # Create output directory
            Path(f"{output_dir}/sketch_samples").mkdir(parents=True, exist_ok=True) 
            Path(f"{output_dir}/gather").mkdir(parents=True, exist_ok=True)        
            if args.start is not None or args.end is not None:
                process_samples_distributed(samples, config, 'sketch', args.threads, args.start, args.end)
                process_samples_distributed(samples, config, 'gather', args.threads, args.start, args.end)
            else:
                # for sample in samples:
                #     run_sketch_for_sample(sample, config)  
                #     run_gather_for_sample(sample, config)
                sketch_success = parallel_process(
                    items=samples,
                    process_func=sketch_wrapper,
                    func_name="sketch",
                    max_workers=args.threads,
                    config=config
                )
                print(f"Sketch: {sketch_success}/{len(samples)} successful")
                gather_success = parallel_process(
                    items=samples,
                    process_func=gather_wrapper,
                    func_name="gather",
                    max_workers=args.threads,
                    config=config
                )
                print(f"Gather: {gather_success}/{len(samples)} successful")


    # These steps should be run only once after all samples are processed
    if args.task in ['all', 'combine']:
        print("\n" + "="*50)
        print("Combining gather results")
        print("="*50)
        config["check_table"] = args.check
        run_gather_all(config)
    
    if args.task in ['all', 'taxonomy']:
        print("\n" + "="*50)
        print("Mapping taxonomy")
        print("="*50)
        run_map_taxonomy(config)
    
    print("\n" + "="*50)
    print("Pipeline execution completed")
    print("="*50)

if __name__ == "__main__":
    main()