import simple_slurm
import argparse
import time
import os

JOBS_DIR = "jobs"  # os.path.join("jobs", "2025-05-14_PnP_GN_TV")
LOGS_DIR = "logs"  # os.path.join("logs", "2025-05-14_PnP_GN_TV")


CMD = """
matlab -nodisplay -nosplash - nodesktop -r 'lamScale=[{lambda_val}]; stepSize={step_size}; run("main_compareContrastSimulatedData.m");'
"""

LAMBDA_VALS = [0.25, 0.5, 1.0, 2.0, 4.0]
STEP_SIZES = [10.0, 1.0, 0.1, 0.01]

KEY_FMT = "lambda_{lambda_val}_step_{step_size}"
JOB_NAME = "2025-05-16_SEAGLE_TV_{key}"


def setup_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test", action="store_true", help="Only submits the first job"
    )
    return parser.parse_args()


def main(args: argparse.Namespace):
    slurm_obj = simple_slurm.Slurm(
        cpus_per_task=16,
        mem="100G",
        partition="general",
        nodelist="g002",
        output=os.path.join(LOGS_DIR, f"{simple_slurm.Slurm.JOB_NAME}.out"),
        error=os.path.join(LOGS_DIR, f"{simple_slurm.Slurm.JOB_NAME}.err"),
    )

    # Create the jobs directory if it doesn't exist
    os.makedirs(JOBS_DIR, exist_ok=True)
    os.makedirs(LOGS_DIR, exist_ok=True)

    # Loop through the lambda values and create job scripts
    for lambda_val in LAMBDA_VALS:
        for step_size in STEP_SIZES:
            slurm_obj.reset_cmd()

            # Get the experiment key
            key = KEY_FMT.format(
                lambda_val=lambda_val,
                step_size=step_size,
            )

            # Set the result directory in the command
            cmd = CMD.format(
                lambda_val=lambda_val,
                step_size=step_size,
            )

            # Reset the slurm object for each job

            # Create the job script
            job_name = JOB_NAME.format(key=key)
            job_path = os.path.join(JOBS_DIR, job_name + ".sh")

            if os.path.exists(job_path):
                print(f"Job script {job_path} already exists. Skipping.")
                continue

            slurm_obj.set_job_name(job_name)
            slurm_obj.set_output(
                os.path.join(LOGS_DIR, f"{simple_slurm.Slurm.JOB_NAME}.out")
            )
            slurm_obj.add_cmd(cmd)

            print(f"Submitting job with job_path={job_path}")
            # print(slurm_obj)
            slurm_obj.sbatch("", job_file=job_path)
            if args.test:
                exit(0)

            time.sleep(1.0)


if __name__ == "__main__":
    args = setup_args()
    main(args)
