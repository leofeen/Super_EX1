from generate_jobs import JOB_CONFIGS

NUMBER_OF_BATCHES = 5
COMMON_BATCH_PATH = r"./batch_"


def compute_and_write_omp_averages() -> None:
    omp_jobs = JOB_CONFIGS['OMP']

    for key, threads in omp_jobs.items():
        m, n = key.split()
        grid_max_time = 0
        for thread in threads:
            time_sum = 0
            iterations_sum = 0
            for i in range(1, NUMBER_OF_BATCHES+1):
                with open(COMMON_BATCH_PATH + str(i) + rf"/result_{m}_{n}_{thread}.txt", "r") as file:
                    time, iterations = map(float, file.readline().split())
                    time_sum += time
                    iterations_sum += iterations

            with open("averages.txt", "a") as file:
                average_time = time_sum/NUMBER_OF_BATCHES
                if grid_max_time == 0:
                    grid_max_time = average_time
                file.write(f"{m} {n} {thread} {average_time:5f} {iterations_sum/NUMBER_OF_BATCHES} {grid_max_time/average_time:5f}\n")


if __name__ == "__main__":
    compute_and_write_omp_averages()