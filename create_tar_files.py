import tarfile
import os
import time
from concurrent.futures import ThreadPoolExecutor
import multiprocessing

base_dir = r"D:\istanbul_wrfouts\23112024"
set_names = [f"SET{i}" for i in range(1, 16)]

def create_tar(set_name):
    source_dir = os.path.join(base_dir, set_name, "temp_extract_dir")
    output_tar = os.path.join(base_dir, set_name, f"{set_name}_wrfouts.tar.gz")

    with tarfile.open(output_tar, "w:gz") as tar:
        for filename in os.listdir(source_dir):
            file_path = os.path.join(source_dir, filename)
            if os.path.isfile(file_path):
                tar.add(file_path, arcname=filename)

    print(f"{set_name}: Done")

# Auto-detect thread count
logical_cores = multiprocessing.cpu_count()
max_workers = min(logical_cores, 8)  # Use up to 8 threads max
print(f"Auto-detected thread count: {max_workers}")

# Start total timing
start_time = time.time()

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    executor.map(create_tar, set_names)

end_time = time.time()
elapsed = end_time - start_time
minutes = int(elapsed // 60)
seconds = int(elapsed % 60)

print(f"All tar.gz archives created with {max_workers} parallel threads.")
print(f"Total elapsed time: {minutes} minutes {seconds} seconds")
