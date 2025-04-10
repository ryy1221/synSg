# This notebook parse gene lists into chopchop to identify all possible sgRNAs
# It will return sgRNA lists for each gene in a folder
import subprocess
import os
import time

# Load gene names from the file
# for name in ['aa','ab','ac','ad','ae','af','ag','ah','ai','aj','ak','al']
for name in ['ae','af','ag','ah','ai','aj','ak','al']:
    file_path = f"../chopchop/whole_genome_genes/mart_export_{name}"
    with open(file_path, "r") as f:
        genes = [line.strip() for line in f if line.strip()]

    print(f"Loaded {len(genes)} genes.")

    # Create directories for scripts and logs
    os.makedirs("slurm_scripts", exist_ok=True)
    os.makedirs("logs", exist_ok=True)

    # Function to check running jobs
    def check_jobs():
        """Check the number of running CHOPCHOP jobs in SLURM."""
        result = subprocess.run(["squeue", "-u", os.getenv("USER")], capture_output=True, text=True)
        return sum(1 for line in result.stdout.split("\n") if "chop_" in line)

    # Function to submit a batch of jobs
    def submit_batch(batch):
        """Submit a batch of jobs using slurmj."""
        for gene in batch:
            job_script = f"""#!/bin/bash
    #SBATCH --job-name=chop_{gene}
    #SBATCH -o /storage/group/epo2/default/yur97/github/synSg/job_status/chopchop_{gene}.out # STDOUT
    #SBATCH --time=1:00:00
    #SBATCH --mem=8G
    #SBATCH --cpus-per-task=1
    #SBATCH --partition=standard

    # Activate Conda environment
    source /swst/apps/anaconda3/2021.05_gcc-8.5.0/etc/profile.d/conda.sh
    conda activate chop

    # Define output paths
    OUTPATH="../data/output_genome/sg_out/{gene}"
    mkdir -p "$OUTPATH"
    OUTF="$OUTPATH/{gene}.txt"

    echo "[`date`] Processing: {gene}"
    ./chopchop.py -p 0 -G hg38 --PAM NGN -o "$OUTPATH" -Target "{gene}" >"$OUTF"
    echo "[`date`] {gene} successful!"
    """

            # Write job script
            script_path = f"slurm_scripts/chopchop_{gene}.slurm"
            with open(script_path, "w") as f:
                f.write(job_script)

            # Submit the job using slurmj
            try:
                subprocess.run(["sbatch", script_path], check=True)
                # print(f"Submitted job for {gene}.")
            except subprocess.CalledProcessError as e:
                print(f"Failed to submit job for {gene}: {e}")

    # Process genes in batches of 80
    batch_size = 60
    for i in range(0, len(genes), batch_size):
        batch = genes[i : i + batch_size]

        print(f"\nSubmitting batch {i // batch_size + 1} ({len(batch)} genes)...")
        submit_batch(batch)

        # Wait for all jobs to finish before submitting the next batch
        print("Waiting for jobs to finish...")
        while check_jobs() > 20:
            print(f"{check_jobs()} jobs still running... checking again in 30 seconds.")
            time.sleep(30)  # Wait 2 minutes before checking again

    print("All jobs completed!")