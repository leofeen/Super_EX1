CORES_PER_NODE = 20

JOB_CONFIGS = {
    "sequential": {
        "20 20": [0],
        "40 40": [0],
        "80 80": [0],
    },

    "OMP": {
        "400 600": [2, 4, 8, 16],
        "800 1200": [4, 8, 16, 32],
    },
}

# SH_RUN_SCRIPT_FILEPATH = r"./run_batch.sh"

def generate_seq_file(m: str, n:str, index: int) -> None:
    with open(rf"./jobs/job_seq_{m}_{n}_{index}.lsf", "w") as file:
        file.write(f'#BSUB -m "polus-c2-ib polus-c4-ib"\n#BSUB -W 00:15\n#BSUB -n 1\n#BSUB -o "SEQ.%J.out"\n#BSUB -e "SEQ.%J.err"\n#BSUB -R "affinity[core(1)]"\n./SEQ {m} {n} {index}')

def generate_omp_file(m: str, n:str, threads: int) -> None:
    with open(rf"./jobs/job_omp_{m}_{n}_{threads}.lsf", "w") as file:
        if threads <= CORES_PER_NODE:
            file.write(f'#BSUB -m "polus-c2-ib polus-c4-ib"\n#BSUB -W 00:15\n#BSUB -n 1\n#BSUB -o "OMP.%J.out"\n#BSUB -e "OMP.%J.err"\n#BSUB -R "affinity[core({threads})]"\nexport OMP_NUM_THREADS={threads} OMP_PROC_BIND=close\n/polusfs/lsf/openmp/launchOpenMP.py ./OMP {m} {n} {threads}')
        else:
            file.write(f'#BSUB -m "polus-c2-ib polus-c4-ib"\n#BSUB -W 00:15\n#BSUB -n 1\n#BSUB -o "OMP.%J.out"\n#BSUB -e "OMP.%J.err"\n#BSUB -R "affinity[core({threads//2})]"\nexport OMP_NUM_THREADS={threads} OMP_PROC_BIND=close\n/polusfs/lsf/openmp/launchOpenMP.py ./OMP {m} {n} {threads}')

def generate_all_files() -> None:
    # with open(SH_RUN_SCRIPT_FILEPATH, "w") as file:
        # file.write("1#!/bin/bash\n")

        for key, _ in JOB_CONFIGS['sequential'].items():
            m, n = key.split()
            generate_seq_file(m, n, 0)
            # file.write(rf'echo "bsub < ./jobs/job_seq_{m}_{n}_0.lsf"' + '\n')

        for key, threads in JOB_CONFIGS['OMP'].items():
            m, n = key.split()
            for thread in threads:
                generate_omp_file(m, n, thread)
                # file.write(rf'echo "bsub < ./jobs/job_omp_{m}_{n}_{thread}.lsf"' + "\n")


if __name__ == "__main__":
    generate_all_files()