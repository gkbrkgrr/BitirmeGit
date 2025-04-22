import os
import tarfile

def listWrfouts(archive_path, extract_dir):
    os.makedirs(extract_dir, exist_ok=True)

    extracted_files = [
        os.path.join(root, file)
        for root, dirs, files in os.walk(extract_dir)
        for file in files if file.startswith("wrfout")
    ]

    if len(extracted_files) != 146:
        with tarfile.open(archive_path, "r:gz") as tar:
            for member in tar.getmembers():
                member.name = member.name.replace(":", "_")  
                tar.extract(member, path=extract_dir)

        extracted_files = [
            os.path.join(root, file)
            for root, dirs, files in os.walk(extract_dir)
            for file in files if file.startswith("wrfout")
        ]

    return sorted(extracted_files)